# ========================================================================== #
# LPJ-LMfire Global Totals Spreadsheet Generator
# ========================================================================== #
# 
# This script generates reproducible spreadsheets of global totals from 
# LPJ-LMfire output data under three main categories of experiments:
# 
#   1. **Master Spreadsheet**:  
#      Global totals from the primary study simulations (fire-only, 
#      photorespiration-only, and both fire & photorespiration).
#
#   2. **Sensitivity Spreadsheets**:  
#      Global totals for parameter uncertainty tests where min/max parameter 
#      combinations are systematically varied.
#
#   3. **Climate Configuration Spreadsheets**:  
#      Global totals under alternative climate scenarios (e.g. high CO2, 
#      high temperature).
#
# Spreadsheets from the study are available in the `data/` directory, but this
# script provides a workflow to regenerate them from the raw LPJ-LMfire output.
#
# ========================================================================== #



# -------------------------------------------------------------------------- #
# Function: generate_global_totals
# -------------------------------------------------------------------------- #
# Purpose:
#   Process LPJ-LMfire output and generate a spreadsheet of global totals 
#   for different oxygen concentrations under multiple simulation scenarios.
#
# Inputs:
#   - basepath    : Directory containing LPJ output files
#   - scriptpath  : Directory containing helper R scripts
#   - simulations : Vector of simulation names (default: fire_only, photo_only, fire_and_photo)
#   - simfolders  : Corresponding subfolders for simulations
#   - config      : Configuration label (for climate configurations)
#   - outfile     : Name of the Excel output file
#   - varnames    : Variables to extract (default: NPP, pftalbiomass, forestcov)
#   - oxygen      : Vector of O2 concentrations to process
#   - save_dir    : Directory to save the Excel output
#
# Outputs:
#   - Excel spreadsheet (.xlsx) of global totals and normalized values.
#
# Notes:
#   - If an oxygen level fails (e.g. missing/corrupt file), it is skipped.
#   - Normalization is relative to O2 = 20.95%.
# -------------------------------------------------------------------------- #

generate_global_totals <- function(
    basepath,
    scriptpath,
    simulations   = c("fire_only", "photo_only", "fire_and_photo"),
    simfolders    = c("oxygen_fire", "oxygen_productivity", "oxygen_fire_productivity"),
    config        = "",
    outfile       = "",
    varnames      = c("NPP", "pftalbiomass", "forestcov"),
    oxygen        = c(16,17,18,19,20,20.95,22,23,24,25,26,27,28,30,31,32,33,34,35),
    save_dir      = getwd()
) {
  # --- Load libraries --- #
  library(ncdf4)
  library(openxlsx)
  
  # --- Load external helper functions --- #
  source(file.path(scriptpath, "decadal_avg.R"))
  source(file.path(scriptpath, "global_total_function.R"))
  
  # --- Set reference O2 (for normalization) --- #
  ref_norm <- which(oxygen == 20.95)
  if (length(ref_norm) != 1) {
    stop("Reference oxygen level 20.95% not found in oxygen vector.")
  }
  
  # --- Initialize results --- #
  result_list <- list()
  valid_oxygen <- oxygen  # track which oxygen levels remain
  
  # --- Loop over simulations --- #
  for (i in seq_along(simulations)) {
    sim <- simulations[i]
    folder <- simfolders[i]
    cat("Processing simulation:", sim, "in folder:", folder, "\n")
    
    tempvar <- list()
    normvals <- list()
    oks_used <- c()  # oxygen values that worked
    
    # --- Loop over variables --- #
    for (v in seq_along(varnames)) {
      var <- varnames[v]
      tempox <- c()
      oks_this_var <- c()
      
      # --- Loop over oxygen levels --- #
      for (o in seq_along(oxygen)) {
        ox <- oxygen[o]
        
        # Get the correct filename for each scenario
        if(sim == "fire_only"){
          ncfile <- file.path(basepath, folder, paste0(ox, "_fire_july22.nc"))
        }else if(sim == "photo_only"){
          ncfile <- file.path(basepath, folder, paste0(ox, "_andre.nc"))
        }else if(sim == "fire_and_photo"){
          ncfile <- file.path(basepath, folder, paste0(ox, "_both_oct.nc"))
        }else{
          ncfile <- file.path(basepath, folder, paste0(config, "_", sim, "_O2", ox, ".nc"))
        }
        
        # Attempt to process file
        val <- tryCatch({
          mdf <- decadalavg_mdf(ncfile, var, path = scriptpath)
          
          if (var == "forestcov") {
            mdf$z.value[mdf$z.value < 0.6] <- 0
            mdf$z.value[mdf$z.value >= 0.6] <- 1
            global_total(mdf, FOREST = TRUE)
          } else {
            global_total(mdf, BIO = TRUE)
          }
        },
        error = function(e) {
          warning(paste("Skipping oxygen level", ox, "for", sim, var, ":", e$message))
          return(NULL)
        })
        
        # Store successful values
        if (!is.null(val)) {
          tempox <- c(tempox, val)
          oks_this_var <- c(oks_this_var, ox)
        }
      }
      
      # Store results for this variable
      tempvar[[v]] <- tempox
      oks_used <- union(oks_used, oks_this_var)
    }
    
    # --- Align variables by common oxygen levels --- #
    oks_final <- Reduce(intersect, lapply(tempvar, function(x) oks_used))
    idx_ref <- which(oks_final == 20.95)
    if (length(idx_ref) == 0) {
      stop("Reference O2=20.95 missing after skipping failed levels!")
    }
    
    mat_vals <- matrix(nrow = length(oks_final), ncol = length(varnames))
    mat_norm <- matrix(nrow = length(oks_final), ncol = length(varnames))
    
    for (v in seq_along(varnames)) {
      ox_filter <- oks_final %in% oks_used
      vals_v <- tempvar[[v]][ox_filter]
      mat_vals[, v] <- vals_v
      mat_norm[, v] <- vals_v / vals_v[idx_ref]
    }
    
    result_list[[sim]] <- list(
      val = mat_vals,
      norm = mat_norm,
      oxygen = oks_final
    )
    
    # Update valid oxygen levels across all simulations
    valid_oxygen <- intersect(valid_oxygen, oks_final)
  }
  
  # --- Build final output dataframe --- #
  df_parts <- list()
  output_order <- c("fire_only", "photo_only", "fire_and_photo")
  
  for (sim in output_order) {
    sim_label <- toupper(sim)
    vals <- result_list[[sim]]$val
    norms <- result_list[[sim]]$norm
    oks   <- result_list[[sim]]$oxygen
    
    keep_idx <- oks %in% valid_oxygen
    vals <- vals[keep_idx, ]
    norms <- norms[keep_idx, ]
    
    cols <- list()
    sim_headers <- c()
    
    for (v in seq_along(varnames)) {
      var <- varnames[v]
      cols[[length(cols) + 1]] <- vals[, v]
      sim_headers <- c(sim_headers, paste0(sim_label, "_", var, "_val"))
      cols[[length(cols) + 1]] <- norms[, v]
      sim_headers <- c(sim_headers, paste0(sim_label, "_", var, "_norm"))
    }
    
    df_part <- as.data.frame(do.call(cbind, cols))
    colnames(df_part) <- sim_headers
    df_parts[[sim]] <- df_part
  }
  
  outdf <- data.frame(O2 = valid_oxygen)
  for (sim in output_order) {
    outdf <- cbind(outdf, df_parts[[sim]])
  }
  
  # --- Save Excel --- #
  outfile_path <- file.path(save_dir, outfile)
  write.xlsx(outdf, outfile_path, rowNames = FALSE)
  
  cat("Done! Output saved as", outfile_path, "\n")
}



# -------------------------------------------------------------------------- #
# Function: generate_param_spreadsheet
# -------------------------------------------------------------------------- #
# Purpose:
#   Generate spreadsheets for parameter sensitivity experiments.
#   Each spreadsheet contains results across oxygen levels and all 
#   parameter combinations (min/max).
#
# Inputs:
#   - simname  : Simulation type ("fire", "photo", "both")
#   - vars     : Variables to extract
#   - scriptpath: Path to helper scripts
#   - basepath : Directory containing parameter test output
#   - O2_concs : Vector of oxygen concentrations to process
#   - save_dir : Directory to save the Excel output
#
# Outputs:
#   - Excel spreadsheet (.xlsx) of totals across parameter combinations.
#
# Notes:
#   - Missing/corrupt files are skipped with warnings.
# -------------------------------------------------------------------------- #

generate_param_spreadsheet <- function(
    simname,
    vars,
    scriptpath,
    basepath,
    O2_concs = c(20.95,23,25,27,29,31,33,35),
    save_dir = getwd()
) {
  # --- Load libraries --- #
  library(openxlsx)
  library(ncdf4)
  
  # --- Load external functions --- #
  source(file.path(scriptpath, "decadal_avg.R"))
  source(file.path(scriptpath, "global_total_function.R"))
  
  # --- Parameter sets by simulation --- #
  if (simname == "fire") {
    params <- c("MoE", "HoC")
  } else if (simname == "photo") {
    params <- c("tau25", "nresp")
  } else {
    params <- c("MoE", "HoC", "tau25", "nresp")
  }
  
  # --- Generate parameter combinations --- #
  generate_combinations <- function(params) {
    param_list <- lapply(params, function(x) c(1, 2))
    names(param_list) <- params
    do.call(expand.grid, param_list)
  }
  
  # --- Storage --- #
  total_dat <- list()
  
  # --- Loop over variables --- #
  for (v in seq_along(vars)) {
    var <- vars[v]
    combinations <- generate_combinations(params)
    col_index_map <- apply(combinations, 1, function(x) paste0(var, "_", paste(x, collapse = "_")))
    
    temp_o2_dat <- matrix(ncol = (nrow(combinations) + 1), nrow = length(O2_concs))
    colnames(temp_o2_dat) <- c("O2_conc", col_index_map)
    
    # --- Loop oxygen levels --- #
    for (o in seq_along(O2_concs)) {
      ox <- O2_concs[o]
      temp_param_dat <- numeric(nrow(combinations))
      
      for (i in 1:nrow(combinations)) {
        combination <- combinations[i, ]
        param_values <- paste(paste0(names(combination), "_", combination), collapse = "_")
        filename <- paste0("May_", simname, "_", ox, "_", param_values, ".nc")
        filepath <- file.path(basepath, filename)
        
        val <- tryCatch({
          mdf <- decadalavg_mdf(filepath, var, path = scriptpath)
          if (var == "forestcov") {
            mdf$z.value[mdf$z.value < 0.6] <- 0
            mdf$z.value[mdf$z.value >= 0.6] <- 1
            global_total(mdf, FOREST = TRUE)
          } else if (var == "mnfire") {
            global_total(mdf, COUNT = TRUE)
          } else {
            global_total(mdf, BIO = TRUE)
          }
        }, error = function(e) {
          warning(paste("Skipping:", filepath, ":", e$message))
          return(NA)
        })
        
        temp_param_dat[i] <- val
      }
      
      temp_o2_dat[o, ] <- c(ox, temp_param_dat)
    }
    
    total_dat[[var]] <- as.data.frame(temp_o2_dat)
  }
  
  # --- Combine all variables --- #
  outdf <- Reduce(function(x, y) cbind(x, y[,-1]), total_dat) # drop duplicate O2 column
  
  # --- Save spreadsheet --- #
  outfile <- paste(simname, "param_generated.xlsx", sep = "_")
  outfile_path <- file.path(save_dir, outfile)
  write.xlsx(outdf, outfile_path, rowNames = FALSE)
  
  cat("Done! Parameter test spreadsheet saved as", outfile_path, "\n")
}



# -------------------------------------------------------------------------- #
# Function: create_totals_spreadsheets
# -------------------------------------------------------------------------- #
# Purpose:
#   Driver function to generate all required spreadsheets (master, sensitivity,
#   and climate). Can be called externally from another script.
#
# Inputs:
#   - main_dir        : Path to the project main directory
#   - data_dir        : Path to LPJ-LMfire data directory
#   - which_spreadsheet: "ALL", "master", "sensitivity", or "climate"
#
# Outputs:
#   - Spreadsheets saved in the appropriate subfolders of data/LPJLMfire_output/totals
# -------------------------------------------------------------------------- #

create_totals_spreadsheets <- function(main_dir, data_dir, which_spreadsheet = "ALL"){
  
  # --- Set paths --- #
  code_dir <- file.path(main_dir, "code/")
  scriptpath <- code_dir
  
  # --- Generate master spreadsheet --- #
  if (which_spreadsheet == "ALL" || which_spreadsheet == "master"){
    generate_global_totals(
      basepath   = data_dir,
      scriptpath = scriptpath,
      outfile    = "master_oxygen_totals_generated.xlsx",
      save_dir   = file.path(main_dir, "data/LPJLMfire_output/totals")
    )
  }
  
  # --- Generate sensitivity spreadsheets --- #
  if (which_spreadsheet == "ALL" || which_spreadsheet == "sensitivity"){
    param_sims <- c("fire", "photo", "both")
    param_vars <- c("pftalbiomass","forestcov")
    
    for (sim in param_sims) {
      generate_param_spreadsheet(
        simname   = sim,
        vars      = param_vars,
        scriptpath = scriptpath,
        basepath   = file.path(data_dir, "May_param_tests", paste0(sim, "_param")),
        save_dir   = file.path(main_dir, "data/LPJLMfire_output/totals/parameter_tests")
      )
    }
  }
  
  # --- Generate climate configuration spreadsheets --- #
  if (which_spreadsheet == "ALL" || which_spreadsheet == "climate"){
    climconfigs <- c("LPJ_high_CO2", "LPJ_high_temp", "LPJ_high_CO2_temp","LPJ_high_CO2_temp_precip")
    
    for ( cc in climconfigs) {
      generate_global_totals(
        basepath   = file.path(data_dir, "late_cret_2025", cc),
        scriptpath = scriptpath,
        simfolders = c("fire_only", "photo_only", "fire_and_photo"),
        config     = cc,
        outfile    = paste(cc,"_generated.xlsx",sep=""),
        save_dir   = file.path(main_dir, "data/LPJLMfire_output/totals/climate_configurations")
      )
    }
  }
}
