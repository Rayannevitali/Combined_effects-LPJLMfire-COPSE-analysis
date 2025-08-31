# ============================================================================== #
# Photosynthesis Module from LPJ (based on BIOME3 model, after Haxeltine & Prentice 1996)
# ------------------------------------------------------------------------------ #
# Function: photomod()
#
# Description:
#   Implements the photosynthesis submodule from LPJ-LMfire, derived from BIOME3.
#   Calculates daily gross photosynthesis (Agd) under given atmospheric CO2 and O2.
#
# Inputs:
#   O2   - Atmospheric oxygen concentration (in % volume, e.g. 21 for 21% O2)
#   CO2  - Atmospheric CO2 concentration as a fractional mole fraction 
#          (e.g. 0.0004 for 400 ppm = 0.04% volume)
#   C4   - Logical (default FALSE). If TRUE, photosynthesis follows the C4 pathway.
#
# Outputs:
#   Agd  - Daily gross photosynthesis [gC/m2/day]
#
# Notes:
#   - All C3 vegetation other than C4 grass is assumed to follow the C3 pathway.
#   - Derived from LPJ-LMfire implementation, consistent with BIOME3 and 
#     Haxeltine & Prentice (1996).
#   - Includes temperature inhibition, quantum efficiency, Rubisco kinetics, 
#     and light limitation.
# ============================================================================== #

photomod <- function(O2, CO2, C4 = FALSE) {
  
  # -------------------------------------------------------------------------- #
  # Model Parameters
  # -------------------------------------------------------------------------- #
  alphaa    = 0.5       # PAR scaling factor (ecosystem vs leaf level)
  alphac3   = 0.08      # Quantum efficiency of CO2 uptake (C3)
  alphac4   = 0.053     # Quantum efficiency of CO2 uptake (C4)
  bc3       = 0.015     # Leaf respiration as fraction of Vmax (C3)
  bc4       = 0.02      # Leaf respiration as fraction of Vmax (C4)
  cmass     = 12.0      # Atomic mass of carbon
  cq        = 4.6e-6    # Conversion factor for radiation (J/m2 → E/m2)
  e0        = 308.56    # Arrhenius parameter
  kc25      = 30.0      # Michaelis constant for CO2 at 25°C
  ko25      = 3e4       # Michaelis constant for O2 at 25°C
  lambdamc3 = 0.8       # Optimal ci:ca ratio (C3)
  lambdamc4 = 0.4       # Optimal ci:ca ratio (C4)
  m         = 25.0      # Parameter in Haxeltine & Prentice Eq. 28
  n0        = 7.15      # Leaf N concentration not involved in photosynthesis
  p         = 1e5       # Atmospheric pressure (Pa)
  po2       = 20.9e3    # O2 partial pressure (Pa)
  q10kc     = 2.1       # Q10 for kc
  q10ko     = 1.2       # Q10 for ko
  q10tau    = 0.57      # Q10 for tau
  t0c3      = 250.0     # Arrhenius base temp for C3 (K)
  t0c4      = 260.0     # Arrhenius base temp for C4 (K)
  tau25     = 2600.0    # Tau at 25°C
  theta     = 0.7       # Co-limitation (shape) parameter
  tk25      = 298.15    # 25°C in Kelvin
  tmc3      = 45.0      # Max temperature for C3 photosynthesis
  tmc4      = 55.0      # Max temperature for C4 photosynthesis
  
  # -------------------------------------------------------------------------- #
  # Environmental Inputs (fixed in this version – can be parameterized later)
  # -------------------------------------------------------------------------- #
  ca      = CO2             # CO2 mole fraction
  dayl    = 11              # Day length (h)
  fpar    = 4.1763796E-03   # Fraction of PAR intercepted
  lambda  = 0.9             # Actual ci:ca ratio
  par     = 7233130.0       # Net photosynthetically active radiation
  temp    = 20              # Temperature (°C)
  
  # Temperature inhibition function parameters
  x1 = 2; x2 = 25; x3 = 30; x4 = 55
  lambdam = 0.9
  
  # -------------------------------------------------------------------------- #
  # Step 1: Calculate absorbed PAR (APAR)
  # -------------------------------------------------------------------------- #
  apar = par * fpar * alphaa   # J/m2/day
  
  # -------------------------------------------------------------------------- #
  # Step 2: Calculate temperature inhibition function (tstress)
  # -------------------------------------------------------------------------- #
  if (temp < x4) {
    k1  = 2. * log((1./0.99) - 1.) / (x1 - x2)
    k2  = (x1 + x2) / 2.
    low = 1. / (1. + exp(k1 * (k2 - temp)))
    
    k3   = log(0.99 / 0.01) / (x4 - x3)
    high = 1. - 0.01 * exp(k3 * (temp - x3))
    
    tstress = low * high
  } else {
    tstress = 0.0
  }
  if (tstress < 1e-2) tstress = 0.0
  
  # -------------------------------------------------------------------------- #
  # Step 3: Rubisco-limited photosynthesis pathway (C3 vs C4 split)
  # -------------------------------------------------------------------------- #
  if (C4) {
    # ---------------------- C4 PHOTOSYNTHESIS ------------------------------- #
    c1 = tstress * alphac4
    if (temp > tmc4) c1 = 0.0
    c2 = 1.0
    b  = bc4
    t0 = t0c4
    
  } else {
    # ---------------------- C3 PHOTOSYNTHESIS ------------------------------- #
    ko  = ko25 * q10ko^((temp - 25.) / 10.)  # Michaelis constant for O2
    kc  = kc25 * q10kc^((temp - 25.) / 10.)  # Michaelis constant for CO2
    
    tau = 0.132 * q10tau^((temp - 25.) / 10.)  # CO2/O2 specificity ratio
    gammastar = ((2/3) * (O2 / tau)) * 0.101325  # CO2 compensation point (Pa)
    
    pa = ca * p      # Ambient CO2 (Pa)
    pi = lambdam * pa  # Intercellular CO2 (Pa)
    
    c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2. * gammastar))
    if (temp > tmc3) c1 = 0.0
    
    c2 = (pi - gammastar) / (pi + kc * (1. + (O2*1e3) / ko))
    b  = bc3
    t0 = t0c3
  }
  
  # -------------------------------------------------------------------------- #
  # Step 4: Calculate Rubisco capacity (Vm)
  # -------------------------------------------------------------------------- #
  s     = (24. / dayl) * b
  sigma = sqrt(max(0., 1. - (c2 - s) / (c2 - theta * s)))
  vm    = (1. / b) * (c1 / c2) * ((2. * theta - 1.) * s - (2. * theta * s - c2) * sigma) * apar * cmass * cq
  
  # -------------------------------------------------------------------------- #
  # Step 5: Calculate PAR-limited and Rubisco-limited rates (JE, JC)
  # -------------------------------------------------------------------------- #
  je = c1 * apar * cmass * cq / dayl  # PAR-limited rate
  jc = c2 * vm / 24.                  # Rubisco-limited rate
  
  if (je < 1e-10 || jc <= 1e-10) {
    agd = 0.0
  } else {
    # Daily gross photosynthesis (Eq. 2, corrected for missing theta)
    agd = (je + jc - sqrt((je + jc)^2. - 4. * theta * je * jc)) / (2. * theta) * dayl
  }
  
  # -------------------------------------------------------------------------- #
  # Step 6: Respiration and Net Photosynthesis (not returned here)
  # -------------------------------------------------------------------------- #
  rd = b * vm
  and = agd - rd
  adt = and + (1. - dayl / 24.) * rd
  
  # -------------------------------------------------------------------------- #
  # Return Gross Photosynthesis
  # -------------------------------------------------------------------------- #
  return(agd)
}
