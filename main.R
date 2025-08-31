# main.R
# -------------------------------------------------------------------------- #
# Script to create the figures from Vitali et al. 2025:
# 'Combined effects of photorespiration and fire strongly regulate atmospheric oxygen levels' 


# --- Plot set up ------------------------------------ # 

# Libraries:
library(ncdf4)
library(ggplot2)
library(gtable)
library(mgcv)
library(reshape2)
library(readxl)
library(dplyr)
library(cowplot)
library(gridExtra)
library(tidyr) 
library(R.matlab)
library(deeptime)
library(Cairo)

# Set Directories (change here)
main_dir <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD/Code/Combined_effects_figures" 
LPJ_data <- "/Volumes/PhD stuff/LPJLMfire/Output/" # where downloaded data is stored. 

# Default directories
code_dir <- paste(main_dir,"/code", sep="")
save_dir <- paste(main_dir,"/Figures", sep="")
data_dir <- paste(main_dir,"/data", sep="")

# Data to download and sort (see comments for each)
Mills_data   <- read_excel(paste(data_dir,"/COPSE_output/proxy_data/Mills_etal_2023_AREPS_O2.xlsx",sep="")) # data can be found here: https://www.annualreviews.org/doi/abs/10.1146/annurev-earth-032320-095425
proxy_dat    <- readMat(paste(data_dir,"/COPSE_output/proxy_data/geochem_data_2020.mat",sep=""))            # data can be downloaded here: https://github.com/bjwmills/SCION/tree/main/data 
Alaska_data  <- read.csv(paste(data_dir,"/LPJLMfire_output/supplementary_data/aggregated_agb_results.csv",sep="")) # data can be downloaded here: https://catalogue.ceda.ac.uk/uuid/af60720c1e404a9e9d2c145d2b2ead4e/ and then aggregated using the aggregated_alaska.R script
HoC_data     <- paste(data_dir,"/LPJLMfire_output/supplementary_data/Compiled_HoC.xlsx",sep="")             # data compiled from Vitali et al.2022, provided here

# Source helper functions
setwd(main_dir)
source(file.path(code_dir, "decadal_avg.R"))
source(file.path(code_dir, "globplot_function.R"))
source(file.path(code_dir, "global_total_function.R"))
source(file.path(code_dir, "plot_helper.R"))
source(file.path(code_dir, "photomod_func.R"))


# --- Common parameters ------------------------------ # 

experiments = c("fire","product","both") 
lab = c("A","B","C","")
lablab = c("(a)", "(b)", "(c)")
oxygen_3 = c(20.95,25,35)


# --- Generate global total spreadsheets ------------- #
# spreadsheets are provided in the data dir, here you can reproduce them with the following file and function. 
# Select which spreadsheets to generate from 'All', 'master', 'sensitivity' or 'climate'
# Note that generating all of the files will take 5-10 minutes.

#source(file.path(code_dir, "global_totals_spreadsheet_generator.R"))
#create_totals_spreadsheets(main_dir, LPJ_data, which_spreadsheet = "ALL")



# --- Main text Figures ------------------------------ # 

# -------------------------------------------------------------------------- #
## Figure 1: Global tree cover from O2 simulations ----
# -------------------------------------------------------------------------- #

fig1 <- plot_oxygen_grid(
  experiments   = c("fire","product","both"),
  folders       = c("oxygen_fire","oxygen_productivity","oxygen_fire_productivity"),
  varname       = "treecov",
  oxygen_levels = c(20.95, 25, 35),
  LPJ_data      = LPJ_data,
  code_dir      = code_dir
)

# Save figure plots:
ggsave(fig1,path = save_dir,filename = "Figure_1.pdf", width = 7.24, height = 4.6, units = c("in"))
ggsave(fig1,path = save_dir,filename = "Figure_1.png", width = 7.24, height = 4.6, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 2: Line plots of global vegetation over O2 from LPJ-LMfire ----
# -------------------------------------------------------------------------- #

# load data:
my_data <- read_excel(paste(data_dir,"/LPJLMfire_output/totals/master_oxygen_totals_original.xlsx",sep=""),sheet = "productivity")
my_data <- data.frame(my_data)

# plot variables:
vplot <- c()
vval  <- c(4,6) # setting which columns in the spreadsheet to select
labs  <- c("biomass","forest cover")
lims  <- c(1400,90)
fig   <- lab[1:2]

for (var in 1:2){
  
  # set title
  if (var==1){
    ytitle <- "aboveground biomass (PgC)"
  }else{
    ytitle <- expression("forest cover (km"^2*" x10"^6*")")
  }
  
  # prediction data for smoothing
  preddat <- data.frame(seq(from = 21, to = 35, by = 0.5))
  colnames(preddat) <- c("Ox")
  
  splines <- matrix(ncol = 3, nrow = length(preddat$Ox))
  ubound <- matrix(ncol = 3, nrow = length(preddat$Ox))
  lbound <- matrix(ncol = 3, nrow = length(preddat$Ox))
  
  for (exp in 1:3){
    
    # load data for this variable and experiment
    exdf = data.frame(my_data$O2[c(7:length(my_data$O2))], as.numeric(my_data[c(7:length(my_data$O2)),(vval[var]+((exp - 1)*6))]))
    colnames(exdf) = c("Ox", "var")
    
    # fit model 
    bsmod <- gam(var ~ s(Ox, k=15, bs = "bs", m = c(3, 2)),
                 data = exdf, method = "REML")
    
    # make predictions and save output
    pred <- predict(bsmod, preddat, se.fit = TRUE, unconditional = FALSE)
    
    # set upper bound 
    ub <- pred$fit +1.96*(pred$se.fit)
    
    # set lower bound
    lb <- pred$fit - 1.96*(pred$se.fit)
    lb[which(lb<0)] = 0
    
    # save data for experiment
    if (var == 2){
      splines[,exp] <- pred$fit/ (10^6)
    }else{
      splines[,exp] <- pred$fit
    }
    ubound[,exp] <- ub
    lbound[,exp] <- lb
    splines <- data.frame(splines)
    ubound <- data.frame(ubound)
    lbound <- data.frame(lbound)
  }
  
  vplot[[var]] <- ggplot()+
    geom_line(aes_(preddat$Ox[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0))],splines[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0)),1],col="fire", linetype = "fire"),linewidth = 0.4)+
    geom_line(aes_(preddat$Ox[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0))],splines[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0)),2],col="photorespiration", linetype = "photorespiration"),linewidth=0.4)+
    geom_line(aes_(preddat$Ox[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0))],splines[c(which(preddat$Ox == 21.0):which(preddat$Ox == 35.0)),3],col="combined", linetype = "combined"),linewidth=0.4)+
    scale_color_manual("",values = c("#3A4454","#F15025","#1A8FE3"))+
    scale_fill_manual("",values = c("#3A4454","#F15025","#1A8FE3"))+
    scale_linetype_manual("",values = c("solid","solid","solid"))+
    geom_hline(yintercept = 0, col = "grey60",linetype = "dotted")+
    ggtitle(paste(fig[var]))+
    scale_x_continuous(limits = c(20.9, 35.1), breaks = seq(21, 35, by = 2), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, lims[var]), breaks = seq(0, lims[var], length.out = 6)) +
    xlab("Oxygen (% of atmosphere)") +
    ylab(ytitle) +
    O2_line_theme 
}

# grab legend from plot
leg <- g_legend(vplot[[2]])

# uncertainty:
# Load and melt biomass data
melted_biomass_p_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/photo_param_master.xlsx",sep=""), 1:5)
melted_biomass_f_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/fire_param_master.xlsx", sep=""), 1:5)
melted_biomass_both_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/both_param_master.xlsx", sep=""), 1:17)

# Compute min, max and loess models for biomass data
photo_loess_models <- compute_min_max_and_loess(melted_biomass_p_data)
fire_loess_models <- compute_min_max_and_loess(melted_biomass_f_data)
both_loess_models <- compute_min_max_and_loess(melted_biomass_both_data)

# Generate a sequence for smoother x-axis values
x_seq <- seq(min(photo_loess_models$min_max_data$O2_conc), max(photo_loess_models$min_max_data$O2_conc), length.out = 500)

# Assuming smoothed data frames for biomass and forest cover are already generated
# Biomass data (no division needed)
smoothed_photo_data <- generate_smoothed_data(photo_loess_models, x_seq)
smoothed_fire_data <- generate_smoothed_data(fire_loess_models, x_seq)
smoothed_both_data <- generate_smoothed_data(both_loess_models, x_seq)

# forest data generation
melted_forest_p_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/photo_param_master.xlsx",sep=""), c(1,7:10))
melted_forest_f_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/fire_param_master.xlsx",sep=""), c(1,7:10))
melted_forest_both_data <- load_and_melt_data(paste(data_dir,"/LPJLMfire_output/totals/parameter_tests/both_param_master.xlsx",sep=""), c(1,19:34))

photo_forest_loess_models <- compute_min_max_and_loess(melted_forest_p_data)
fire_forest_loess_models <- compute_min_max_and_loess(melted_forest_f_data)
both_forest_loess_models <- compute_min_max_and_loess(melted_forest_both_data)

# smooth forest data
smoothed_forest_photo_data <- generate_smoothed_data(photo_forest_loess_models, x_seq)
smoothed_forest_fire_data <- generate_smoothed_data(fire_forest_loess_models, x_seq)
smoothed_forest_both_data <- generate_smoothed_data(both_forest_loess_models, x_seq)

# Forest cover data (divide by 1,000,000)
smoothed_forest_photo_data_adjusted <- smoothed_forest_photo_data
smoothed_forest_photo_data_adjusted$min_smooth <- smoothed_forest_photo_data$min_smooth / 1000000
smoothed_forest_photo_data_adjusted$max_smooth <- smoothed_forest_photo_data$max_smooth / 1000000

smoothed_forest_fire_data_adjusted <- smoothed_forest_fire_data
smoothed_forest_fire_data_adjusted$min_smooth <- smoothed_forest_fire_data$min_smooth / 1000000
smoothed_forest_fire_data_adjusted$max_smooth <- smoothed_forest_fire_data$max_smooth / 1000000

smoothed_forest_both_data_adjusted <- smoothed_forest_both_data
smoothed_forest_both_data_adjusted$min_smooth <- smoothed_forest_both_data$min_smooth / 1000000
smoothed_forest_both_data_adjusted$max_smooth <- smoothed_forest_both_data$max_smooth / 1000000

# Create the ribbon layers for biomass data
smooth_ribbon_layer <- geom_ribbon(data = smoothed_photo_data, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#1A8FE3", alpha = 0.09)
smooth_ribbon_f_layer <- geom_ribbon(data = smoothed_fire_data, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#F15025", alpha = 0.09)
smooth_ribbon_b_layer <- geom_ribbon(data = smoothed_both_data, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#3A4454", alpha = 0.07)

# Create the ribbon layers for forest cover data
smooth_forest_ribbon_layer <- geom_ribbon(data = smoothed_forest_photo_data_adjusted, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#1A8FE3", alpha = 0.09)
smooth_forest_ribbon_f_layer <- geom_ribbon(data = smoothed_forest_fire_data_adjusted, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#F15025", alpha = 0.09)
smooth_forest_ribbon_b_layer <- geom_ribbon(data = smoothed_forest_both_data_adjusted, aes(x = O2_conc, ymin = min_smooth, ymax = max_smooth), fill = "#3A4454", alpha = 0.07)

# Sensitivity lines on plot here instead of ribbons
biomass_p_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/photo_param_master.xlsx", sep=""))[1:5]
biomass_f_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/fire_param_master.xlsx", sep=""))[1:5]
biomass_both_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/both_param_master.xlsx", sep=""))[1:17]

forest_p_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/photo_param_master.xlsx", sep=""))[c(1,7:10)]
forest_f_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/fire_param_master.xlsx", sep=""))[c(1,7:10)]
forest_both_data <- read_excel(paste(data_dir, "/LPJLMfire_output/totals/parameter_tests/both_param_master.xlsx", sep=""))[c(1,19:34)]

# --- Prepare biomass data ---
long_biomass_combined <- biomass_both_data %>%
  pivot_longer(cols = starts_with("pftalbiomass"),
               names_to = "biomass_type",
               values_to = "biomass_value") %>%
  mutate(group = "combined")

long_biomass_photo <- biomass_p_data %>%
  pivot_longer(cols = starts_with("pftalbiomass"),
               names_to = "biomass_type",
               values_to = "biomass_value") %>%
  mutate(group = "photorespiration")

long_biomass_fire <- biomass_f_data %>%
  pivot_longer(cols = starts_with("pftalbiomass"),
               names_to = "biomass_type",
               values_to = "biomass_value") %>%
  mutate(group = "fire")

all_biomass_long <- bind_rows(long_biomass_combined, long_biomass_photo, long_biomass_fire)

all_biomass_long <- all_biomass_long %>%
  mutate(alpha = case_when(
    group == "combined" ~ 0.15,
    TRUE ~ 0.4
  ))

# --- Plot biomass ---
final_plot_biomass <- vplot[[1]] +
  geom_line(data = all_biomass_long,
            aes(x = O2_conc, y = biomass_value, group = interaction(group, biomass_type),
                color = group, alpha = alpha),
            linewdith = 0.2, linetype = "dashed") +
  smooth_ribbon_b_layer + smooth_ribbon_layer + smooth_ribbon_f_layer +
  scale_alpha_identity()


# --- Prepare forest cover data ---
long_forest_combined <- forest_both_data %>%
  pivot_longer(cols = starts_with("forestcov"),
               names_to = "forestcov_type",
               values_to = "forestcov_value") %>%
  mutate(forestcov_value = forestcov_value / 1e6,
         group = "combined")

long_forest_photo <- forest_p_data %>%
  pivot_longer(cols = starts_with("forestcov"),
               names_to = "forestcov_type",
               values_to = "forestcov_value") %>%
  mutate(forestcov_value = forestcov_value / 1e6,
         group = "photorespiration")

long_forest_fire <- forest_f_data %>%
  pivot_longer(cols = starts_with("forestcov"),
               names_to = "forestcov_type",
               values_to = "forestcov_value") %>%
  mutate(forestcov_value = forestcov_value / 1e6,
         group = "fire")

all_forest_long <- bind_rows(long_forest_combined, long_forest_photo, long_forest_fire)

all_forest_long <- all_forest_long %>%
  mutate(alpha = case_when(
    group == "combined" ~ 0.15,
    TRUE ~ 0.4
  ))

# --- Plot forest cover ---
final_plot_forest <- vplot[[2]] +
  geom_line(data = all_forest_long,
            aes(x = O2_conc, y = forestcov_value, group = interaction(group, forestcov_type),
                color = group, alpha = alpha),
            linewidth = 0.2, linetype = "dashed") +
  smooth_forest_ribbon_b_layer + smooth_forest_ribbon_layer + smooth_forest_ribbon_f_layer +
  scale_alpha_identity()


# Combine plots:
fig2 =  plot_grid(final_plot_biomass+ theme(legend.position = "none",plot.margin=unit(c(0.1,0.15,0,0.2),"cm")),final_plot_forest+ theme(legend.position = "none",plot.margin=unit(c(0.1,0.1,0,0),"cm")),
                    rel_heights = c(1,1,0.2),leg,axis = 'l', ncol = 1, align = "vh")

# Save figure plots:
ggsave(fig2,path = save_dir,filename = "Figure_2.pdf", width = 3.55, height = 5.3, units = c("in"))
ggsave(fig2,path = save_dir,filename = "Figure_2.png", width = 3.55, height = 5.3, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 3: Global Tree cover under high CO2, Temp and precipitation ----
# -------------------------------------------------------------------------- #

fig3 <- plot_oxygen_grid(
  experiments = c("fire","product","both"),
  folders     = c("fire_only","photo_only","fire_and_photo"),
  varname     = "treecov",
  oxygen_levels = c(20.95, 25, 35),
  LPJ_data    = paste0(LPJ_data,"late_cret_2025/"),
  code_dir    = code_dir,
  climate_config = "LPJ_high_CO2_temp_precip"
)


# Save figure plots:
ggsave(fig3,path = save_dir,filename = "Figure_3.pdf", width = 7.24, height = 4.6, units = c("in"))
ggsave(fig3,path = save_dir,filename = "Figure_3.png", width = 7.24, height = 4.6, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 4:  Line plots of global vegetation under high climate ----
# -------------------------------------------------------------------------- #

fig4 <- make_climate_oxline_plots(
  climate_config = "LPJ_high_CO2_temp_precip",
  short_name = "high_CO2_temp_precip",
  O2_line_theme = O2_line_theme,
  base_dir = data_dir
)

# Build plot grid in main script
main_clim_oxline_plot <- cowplot::plot_grid(
  fig4$vplot[[1]] + theme(legend.position = "none"),
  fig4$vplot[[2]] + theme(legend.position = "none"),
  nrow = 2, align = "hv"
)

# Add legend
fig4 <- cowplot::plot_grid(main_clim_oxline_plot, fig4$legend,
                                            ncol = 1, rel_heights = c(1, 0.15))
print(fig4)

# Save figure plots:
ggsave(fig4,path = save_dir,filename = "Figure_4.pdf", width = 3.55, height = 5.3, units = c("in"))
ggsave(fig4,path = save_dir,filename = "Figure_4.png", width = 3.55, height = 5.3, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 5: COPSE O2, VEG & CO2 over time ----
# -------------------------------------------------------------------------- #

# --- Load data --- #
bg_data_path <- paste(data_dir,"/COPSE_output/",sep="")

setwd(bg_data_path)
none     <- read.csv(paste(bg_data_path,"output/COPSE_no_feedbacks.txt", sep=""))
both     <- read.csv(paste(bg_data_path,"output/COPSE_both_feedbacks.txt", sep=""))
combined <- read.csv(paste(bg_data_path,"output/COPSE_combined_feedbacks.txt", sep=""))
unf      <- read.csv(paste(bg_data_path,"output/COPSE_updated_no_fire_feedbacks.txt", sep=""))
unp      <- read.csv(paste(bg_data_path,"output/COPSE_updated_no_prod_feedbacks.txt", sep=""))

mills_o2 <- Mills_data

# --- COPSE O2 plot --- #
COPSEO2 <- ggplot() +
  geom_ribbon(data = mills_o2, aes(x = Time*-1, ymin = Min, ymax = Max), fill = "grey90", alpha = 0.5) +
  geom_hline(yintercept = 20.95, linetype="dotdash", color="grey70", size=0.4) +
  geom_line(aes(both$time_myr*-1, both$mrO2*100, color = "Fire & \nphotorespiration", linetype = "Fire & \nphotorespiration"), size=0.4) +
  geom_line(aes(none$time_myr*-1, none$mrO2*100, color = "No vegetation-based \noxygen feedback", linetype = "No vegetation-based \noxygen feedback"), size=0.4) +
  geom_line(aes(unf$time_myr*-1, unf$mrO2*100, color="Photorespiration only \n(no fire)", linetype = "Photorespiration only \n(no fire)"), size=0.4) +
  geom_line(aes(unp$time_myr*-1, unp$mrO2*100, color="Fire only \n(no photorespiration)", linetype = "Fire only \n(no photorespiration)"), size=0.4) +
  geom_line(aes(combined$time_myr*-1, combined$mrO2*100, color = "Combined", linetype = "Combined"), size=0.4) +
  scale_linetype_manual(name="feedbacks included:", values=copse_linetypes) +
  scale_color_manual(name ="feedbacks included:", values=copse_colors) +
  scale_x_reverse() +
  xlab("Time (myr)") +
  ylab(expression("O"[2]~"(%)")) +
  coord_geo(dat = list("eras", "periods"), xlim = c(400, 0), ylim=c(10,45),
            pos = list("b", "b"), size = list(2.5,2.5), abbrv = list(FALSE, TRUE), height=list(unit(0.7, "lines"))) +
  copse_theme +
  theme(legend.position="right", legend.title=element_text(size=7))

# --- COPSE VEG plot --- #
COPSEVEG <- ggplot() +
  geom_hline(yintercept = 1, linetype="dotdash", color="grey70", size=0.4) +
  geom_line(aes(both$time_myr*-1, both$VEG, color = "Fire & \nphotorespiration", linetype = "Fire & \nphotorespiration"), size=0.4) +
  geom_line(aes(none$time_myr*-1, none$VEG, color = "No vegetation-based \noxygen feedback", linetype = "No vegetation-based \noxygen feedback"), size=0.4) +
  geom_line(aes(unf$time_myr*-1, unf$VEG, color="Photorespiration only \n(no fire)", linetype = "Photorespiration only \n(no fire)"), size=0.4) +
  geom_line(aes(unp$time_myr*-1, unp$VEG, color="Fire only \n(no photorespiration)", linetype = "Fire only \n(no photorespiration)"), size=0.4) +
  geom_line(aes(combined$time_myr*-1, combined$VEG, color = "Combined", linetype = "Combined"), size=0.4) +
  scale_linetype_manual(name="feedbacks included:", values=copse_linetypes) +
  scale_color_manual(name ="feedbacks included:", values=copse_colors) +
  scale_x_reverse() +
  xlab("Time (myr)") +
  ylab("VEG (normalised)") +
  coord_geo(dat = list("eras", "periods"), xlim = c(400,0), ylim=c(0,2),
            pos=list("b","b"), size=list(2.5,2.5), abbrv=list(FALSE,TRUE), height=list(unit(0.7,"lines"))) +
  copse_theme +
  theme(legend.position="bottom", legend.title.align=0.5, legend.title.position="top")

# --- COPSE CO2 plot --- #
COPSECO2 <- ggplot() +
  geom_hline(yintercept = 280, linetype="dotdash", color="grey70", size=0.4) +
  geom_point(aes(proxy_dat$alkenone.age*-1, proxy_dat$alkenone.co2, color="alkenone", shape="alkenone"), size=2, col="grey90", shape=15) +
  geom_point(aes(proxy_dat$boron.age*-1, proxy_dat$boron.co2, color="boron", shape="boron"), size=2, col="grey90", shape=15) +
  geom_point(aes(proxy_dat$liverwort.age*-1, proxy_dat$liverwort.co2, color="liverwort", shape="liverwort"), size=2, col="grey90", shape=15) +
  geom_point(aes(proxy_dat$paleosol.age*-1, proxy_dat$paleosol.co2, color="paleosol", shape="paleosol"), size=2, col="grey90", shape=15) +
  geom_point(aes(proxy_dat$phytane.age*-1, proxy_dat$phytane.co2, color="phytane", shape="phytane"), size=2, col="grey90", shape=15) +
  geom_point(aes(proxy_dat$stomata.age*-1, proxy_dat$stomata.co2, color="stomata", shape="stomata"), size=2, col="grey90", shape=15) +
  geom_line(aes(both$time_myr*-1, both$RCO2*280, color = "Fire & \nphotorespiration", linetype = "Fire & \nphotorespiration"), size=0.4) +
  geom_line(aes(none$time_myr*-1, none$RCO2*280, color = "No vegetation-based \noxygen feedback", linetype = "No vegetation-based \noxygen feedback"), size=0.4) +
  geom_line(aes(unf$time_myr*-1, unf$RCO2*280, color="Photorespiration only \n(no fire)", linetype = "Photorespiration only \n(no fire)"), size=0.4) +
  geom_line(aes(unp$time_myr*-1, unp$RCO2*280, color="Fire only \n(no photorespiration)", linetype = "Fire only \n(no photorespiration)"), size=0.4) +
  geom_line(aes(combined$time_myr*-1, combined$RCO2*280, color = "Combined", linetype = "Combined"), size=0.4) +
  scale_linetype_manual(name="feedbacks included:", values=copse_linetypes) +
  scale_color_manual(name ="feedbacks included:", values=copse_colors) +
  scale_x_reverse() +
  xlab("Time (Myr)") +
  ylab(expression("CO"[2]~"(ppm)")) +
  coord_geo(dat = list("eras", "periods"), xlim = c(400,0), ylim=c(0,4000),
            pos=list("b","b"), size=list(2.5,2.5), abbrv=list(FALSE,TRUE), height=list(unit(0.7,"lines"))) +
  copse_theme +
  theme(legend.position="bottom", legend.title.align=0.5, legend.title.position="top")

# --- Combine plots --- #
copseleg <- g_legend(COPSEVEG)

copseplots <- plot_grid(
  COPSEO2 + theme(legend.position='none'),
  COPSEVEG + theme(legend.position='none'),
  COPSECO2 + theme(legend.position='none'),
  ncol = 1, align = "v", labels = c("A","B","C"), label_size = 10
)

fig5 <- plot_grid(copseplots, copseleg, ncol = 1, rel_heights = c(10,1))
print(fig5)

# Save figure plots:
ggsave(fig5,path = save_dir,filename = "Figure_5.pdf", width = 6.5, height = 8, units = c("in"))
ggsave(fig5,path = save_dir,filename = "Figure_5.png", width = 6.5, height = 8, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 6: Compensation point and gross photosynthesis in LPJLM-fire ----
# -------------------------------------------------------------------------- #

# --- helper functions for compensation points (Nisbet et al.) ---
ccomp <- function(O2) 2.13*O2 + 3.89
ocomp <- function(CO2) 0.0246*CO2*1e6 + 17.94
comp  <- function(O2)  ((2/3) * (O2 / (0.132*0.57^((20 - 25)/10)))) * 1.01325

# --- parameter grid ---
O2  <- seq(1, 100, 1)
CO2 <- seq(1e-5, 1e-3, 1e-5)   # skip zero
grid <- expand.grid(O2 = O2, CO2 = CO2)

# --- calculate photosynthesis for each O2/CO2 combination ---
grid$photo <- mapply(photomod, grid$O2, grid$CO2)

# --- classify response values into breaks ---
grid$breaks <- cut(grid$photo, breaks = c(-Inf, -1e-5, 1e-5, 0.01, 0.02, Inf))

# --- compensation point lines ---
cline <- data.frame(O2 = O2,  CO2 = ccomp(O2))
sline <- data.frame(O2 = O2,  CO2 = comp(O2))
oline <- data.frame(CO2 = CO2, O2 = ocomp(CO2))

# --- plot ---
fig6 <- ggplot(grid, aes(x = O2, y = CO2, fill = photo)) +
  geom_raster(interpolate = TRUE) +
  geom_line(data = sline, aes(x = O2, y = CO2/(0.9*1e6)), # compensation line (gamma*)
            linewidth = 0.4, inherit.aes = FALSE) +
  annotate("text", x = 95, y = 0.00048, label = expression(gamma^"*"), size = 5) +
  geom_vline(xintercept = 20.95, linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = 2.8813e-4, linetype = "dotted", linewidth = 0.3) +
  scale_fill_gradient2(
    name = expression("gross photosynthesis (gC m"^-2*" day"^-1*")"),
    low = "red", high = "green3",
    guide = guide_colourbar(direction = "horizontal", 
                            title.position = "top", title.hjust = 0.5)) +
  scale_y_continuous(
    breaks = seq(0, 1e-3, 1e-4),
    labels = seq(0, 1000, 100)) +
  labs(x = expression(O[2]*" (%)"), y = expression(CO[2]*" (ppm)")) +
  theme(
    panel.border = element_rect(fill=NA, colour="black", linewidth=0.1),
    panel.grid     = element_blank(),
    panel.background = element_blank(),
    axis.text      = element_text(size = 8),
    axis.title     = element_text(size = 8),
    axis.title.y   = element_text(margin = margin(r = 5)),
    axis.title.x   = element_text(margin = margin(t = 5)),
    legend.title   = element_text(size = 7),
    legend.text    = element_text(size = 7),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.spacing.x = unit(0.75, 'cm'),
    legend.margin   = margin(5,0,0,0),
    legend.box.margin = margin(-10,-10,0,-10),
    plot.margin   = unit(c(0.1,0.1,0.1,0.1), "cm")
  )

print(fig6)

# Save figure plots:
ggsave(fig6,path = save_dir,filename = "Figure_6.pdf", width = 3.45, height = 4, units = c("in"))
ggsave(fig6,path = save_dir,filename = "Figure_6.png", width = 3.45, height = 4, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure 7: Feedback Mechanisms and Normalised biomass ----
# -------------------------------------------------------------------------- #

# Plotting panel B:
my_data <- read_excel(paste(data_dir,"/LPJLMfire_output/totals/master_oxygen_totals_original.xlsx",sep="")) # load spreadsheet data 

# interpolation:
fbion <- smooth.spline(as.numeric(unlist(my_data[3:21,1])), 
                       as.numeric(unlist(my_data[3:21,5])),
                       spar = 0.5)

pbion <- smooth.spline(as.numeric(unlist(my_data[3:21,1])), 
                       as.numeric(unlist(my_data[3:21,11])), 
                       spar = 0.5)

bbion <- smooth.spline(as.numeric(unlist(my_data[3:21,1])), 
                       as.numeric(unlist(my_data[3:21,17])), 
                       spar = 0.5)

# plot normalised values
fig7 <- ggplot() +
  geom_line(aes(fbion$x, fbion$y, color = "fire"), size = 0.4) +
  geom_line(aes(pbion$x, pbion$y, color = "photo"), size = 0.4) +
  geom_line(aes(bbion$x, bbion$y, color = "combined"), size = 0.4) +
  scale_y_continuous(breaks = seq(0, 2.7, by = 0.4), limits = c(0, 1.8)) +
  scale_x_continuous(limits = c(15, 35)) +
  scale_color_manual(
    name = "oxygen effects:",
    values = c(
      "fire"     = "#F15025",  
      "photo"    = "#1A8FE3",  
      "combined" = "#3A4454"   
    ),
    labels = c(
      fire     = expression("fire " * italic("(") * italic(f)[italic(fireLPJ)] * italic(")")),
      photo    = expression("photo " * italic("(") * italic(V)[O2,LPJ] * italic(")")),
      combined = expression("combined " * italic("(") * italic(V)[combined] * italic(")"))
    )
  ) +
  xlab("Atmospheric oxygen concentration (%)") +
  ylab("terrestrial biomass (normalised)") +
  O2_line_theme

print(fig7)

# Save figure plots:
ggsave(fig7,path = save_dir,filename = "Figure_7b.pdf", width = 5.5, height = 3, units = c("in"))
ggsave(fig7,path = save_dir,filename = "Figure_7b.png", width = 5.5, height = 3, units = c("in"),dpi=900)




# --- Supplementary Figures --------------------------- # 

# -------------------------------------------------------------------------- #
## Figure S1: Negative feedback diagram ----
# -------------------------------------------------------------------------- #

# Note that this figure was manually created

# -------------------------------------------------------------------------- #
## Figure S2: Sensitivity parameter test range ----
# -------------------------------------------------------------------------- #

# Moisture of extinction function from Watson & Lovelock
MoE <- function(Ox) {
  (8 * Ox) - 128
}

# Plot linear bounds
plot_bounds <- function(x_vals, y1, y2, col='red') {
  upper_slope <- (max(y2, na.rm=TRUE) - max(y1, na.rm=TRUE)) / (max(x_vals) - min(x_vals))
  lower_slope <- (min(y2, na.rm=TRUE) - min(y1, na.rm=TRUE)) / (max(x_vals) - min(x_vals))
  
  lines(c(min(x_vals), max(x_vals)), c(max(y1, na.rm=TRUE), max(y2, na.rm=TRUE)), 
        type='l', lty='dashed', col=col)
  lines(c(min(x_vals), max(x_vals)), c(min(y1, na.rm=TRUE), min(y2, na.rm=TRUE)), 
        type='l', lty='dashed', col=col)
}

# Load data
data_file <- paste(data_dir,"/LPJLMfire_output/supplementary_data/Compiled_HoC.xlsx",sep="")
compiled_data <- read_excel(data_file)

amb <- compiled_data$`Ambient HoC (kJ/g)`
bomb <- compiled_data$`Bomb HoC (kJ/g)`

# Set up 2-panel layout with equal spacing
par(mfrow=c(2,1),
    mar=c(3.5,3.5,1.5,1.5),   # bottom, left, top, right margins per plot
    oma=c(2,2,2,2),       # outer margins for outer labels
    mgp=c(2,0.5,0))       # axis title, labels, line positions


axis_text_size <- 0.8  # roughly matches ggplot2 size=8
axis_label_size <- 0.8
tag_cex <- 1.2

# Panel A: MoE
Ox <- 21:35
M <- MoE(Ox)

plot(Ox, M, type='l', xlim=c(10,40), ylim=c(0,170),
     xlab="Oxygen (% of the atmosphere)", 
     ylab="MoE (%)",
     cex.axis=axis_text_size, 
     cex.lab=axis_label_size, 
     lwd=1,
     font.lab=1)

# Add points and dashed lines
ex <- c(12, 18, 35, 35)
ey <- c(0, 0, 60, 155)
points(ex, ey)
lines(c(12,35), c(0,155), lty='dashed', col='red', lwd=1)
lines(c(18,35), c(0,60), lty='dashed', col='red', lwd=1)

# Label
mtext("A", side=3, line=-1, adj=-0.1, outer=TRUE, cex=tag_cex, font=2)

# Panel B: Heat of Combustion
plot(c(rep(21, length(amb)), rep(100, length(bomb))),
     c(amb, bomb),
     xlab="Oxygen (% of the atmosphere)",
     ylab="Heat of Combustion (kJ/g)",
     main="",
     cex.axis=axis_text_size, 
     cex.lab=axis_label_size)

plot_bounds(c(21,100), amb, bomb)

# Label
mtext("B", side=3, line=-17, adj=-0.1, outer=TRUE, cex=tag_cex, font=2)

# Reset plotting
par(mfrow=c(1,1))


# -------------------------------------------------------------------------- #
## Figure S3: Tree height over atmospheric oxygen concentrations ----
# -------------------------------------------------------------------------- #

# set parameters and arrays
oxygen = c(20.95,35)
var = "height"
heightmdf = c()
covmdf = c()
hfrac = c()
gplot=c()
gp=1

# loop over experiments
for (exp in 1:length(experiments)) {
  ex = experiments[exp]
  
  if( ex == "fire"){
    fname = "_fire_july22.nc"
    pname = "oxygen_fire"
  }else if(ex == "product"){
    fname="_andre.nc"
    pname = "oxygen_productivity"
    ex= "photorespiration"
  }else{
    fname="_both_oct.nc"
    pname = "oxygen_fire_productivity"
    ex= "combined"}
  
  # set temporary array to hold global totals if needed
  tempox = array(0,length(oxygen))
  
  # loop over oxygen concentration
  for (o in 1:length(oxygen)) {
    ox = oxygen[o]
    
    # set correct working directory
    setwd(paste(LPJ_data,pname,sep=""))
    
    # get correct file
    oxfile=(paste(ox,fname,sep=""))
    
    for (pft in 1:7){
      covmdf[[pft]] = decadalavg_mdf(oxfile,"cover",npft = pft, path = code_dir)
      heightmdf[[pft]] = decadalavg_mdf(oxfile,"height",npft=pft, path = code_dir)
      
      # calculate total tree cover in each cell (overall)
      if(pft == 1){
        totcovmdf = covmdf[[1]]
      }else{
        totcovmdf = totcovmdf + covmdf[[pft]]
      }
    }
    
    for(pft in 1:7){
      # calculate fraction of each tree PFT cover out of total tree cover
      covmdf[[pft]] = covmdf[[pft]]/totcovmdf
      covmdf[[pft]][which(is.na(covmdf[[pft]][,3])),3]=0
      
      # multiple fraction against average height for each pft
      hfrac[[pft]] = heightmdf[[pft]] * covmdf[[pft]]
      
      # sum across all pfts to get gridcell average tree height using weighted fractions 
      if(pft==1){
        hmdf = hfrac[[1]]
      }else{
        hmdf = hmdf + hfrac[[pft]]
      }
    }
    
    #create global plot
    labs = paste(lab[exp],lablab[o], sep="")
    gplot[[gp]] = globalplot(hmdf,diff=FALSE,sim = paste(ex),plotlab=labs,o2=ox, var=var)
    
    gp= gp+1 # increase counter
  } 
}

# grab legend from plot
leg = gtable_filter(ggplot_gtable(ggplot_build(gplot[[1]])), "guide-box")

figS3 = grid.arrange(
  arrangeGrob(gplot[[1]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),gplot[[2]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),
              gplot[[3]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),gplot[[4]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),
              gplot[[5]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),gplot[[6]]+ theme(legend.position = "none",plot.margin=unit(c(-0.5,0.1,-1,0),"cm")),
              nrow = 3),leg,heights = c(10,0.7),top="")

# Save figure plots:
ggsave(figS3,path = save_dir,filename = "Figure_S3.pdf", width = 7.24, height = 6.2, units = c("in"))
ggsave(figS3,path = save_dir,filename = "Figure_S3.png", width = 7.24, height = 6.2, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure S4: COPSE model comparison to the geological record ----
# -------------------------------------------------------------------------- #
setwd(bg_data_path)

# load COPSE output
COPSE_stand <- read.csv(paste(bg_data_path,"output/COPSE_standard.txt", sep=""))
COPSE_comb  <- combined

# CO2 proxy data list
CO2_proxy = c("alkenone.co2","boron.co2","liverwort.co2","paleosol.co2","phytane.co2","stomata.co2")

# CO2 plot with coloured proxies
CO2 = ggplot()+
  geom_point(aes(proxy_dat$alkenone.age*-1,proxy_dat$alkenone.co2,color="alkenone", shape="alkenone"),size=2)+
  geom_point(aes(proxy_dat$boron.age*-1,proxy_dat$boron.co2,color="boron", shape="boron"),size=2)+
  geom_point(aes(proxy_dat$liverwort.age*-1,proxy_dat$liverwort.co2,color="liverwort", shape="liverwort"),size=2)+
  geom_point(aes(proxy_dat$paleosol.age*-1,proxy_dat$paleosol.co2,color="paleosol", shape="paleosol"),size=2)+
  geom_point(aes(proxy_dat$phytane.age*-1,proxy_dat$phytane.co2,color="phytane", shape="phytane"),size=2)+
  geom_point(aes(proxy_dat$stomata.age*-1,proxy_dat$stomata.co2,color="stomata", shape="stomata"),size=2)+
  geom_line(aes(COPSE_comb$time_myr*-1, COPSE_comb$RCO2*280))+
  geom_line(aes(COPSE_stand$time_myr*-1, COPSE_stand$RCO2*280), lty='dashed')+ #PAL PPM CO2
  scale_x_reverse()+
  labs(
    x = "Time (Myr)", 
    y = expression("CO"[2]~"(ppm)"),  # Use expression for proper subscripts
    title = "A"
  ) +
  coord_geo(dat = list("periods", "eras"),xlim = c(450, 0),ylim=c(0,4000),pos = list("b", "b"),
            height = unit(0.7, "line"), size = "auto", abbrv = list(TRUE, FALSE)) +
  scale_color_manual("Proxy",values = c(
    "alkenone" = "#E69F00",  # Bright yellow-orange
    "boron" = "#56B4E9",     # Bright sky blue
    "liverwort" = "#009E73", # Bright green
    "paleosol" = "#F0E442",  # Bright yellow
    "phytane" = "#D55E00",   # Bright red-orange
    "stomata" = "#CC79A7"    # Bright pink
  )) + 
  scale_shape_manual("Proxy",values = c(
    "alkenone" = 16, # Circle
    "boron" = 17,    # Triangle
    "liverwort" = 18, # Square
    "paleosol" = 15,  # Diamond
    "phytane" = 19,   # Circle with border
    "stomata" = 20    # Filled circle
  )) +# Custom shapes
  guides(color = guide_legend(ncol = 6)) + # Arrange legend items in 1 row (adjust 'ncol' as needed)
  proxy_theme

# plot 
D13C= ggplot()+
  geom_point(aes(proxy_dat$d13c.x[1,]*-1,proxy_dat$d13c.y[1,]),size=2, col="gold")+
  geom_line(aes(COPSE_comb$time_myr*-1, COPSE_comb$d13c_A))+
  geom_line(aes(COPSE_stand$time_myr*-1, COPSE_stand$d13c_A), lty='dashed')+ #PAL PPM CO2
  scale_x_reverse()+
  labs(
    x = "Time (Myr)", 
    y = expression(paste(delta^{13}, "C"["carb"], " ", "‰")),  # Use "‰" for per mille
    title = "B"
  ) +
  coord_geo(dat = list("periods", "eras"),xlim = c(450, 0),ylim=c(-8,10),pos = list("b", "b"),
            height = unit(0.7, "line"), size = "auto", abbrv = list(TRUE, FALSE)) +
  guides(color = guide_legend(ncol = 6)) + # Arrange legend items in 1 row (adjust 'ncol' as needed)
  proxy_theme


D34S= ggplot()+
  geom_point(aes(proxy_dat$d34s.x[1,]*-1,proxy_dat$d34s.y[1,]),size=2, col="gold")+
  geom_line(aes(COPSE_comb$time_myr*-1, COPSE_comb$d34s_S))+
  geom_line(aes(COPSE_stand$time_myr*-1, COPSE_stand$d34s_S), lty='dashed')+ #PAL PPM CO2
  scale_x_reverse()+
  labs(
    x = "Time (Myr)", 
    y = expression(paste(delta^{34}, "S"["SW"], " ", "‰")),  # Use "‰" for per mille
    title = "C"
  ) +
  coord_geo(dat = list("periods", "eras"),xlim = c(450, 0),ylim=c(-10,50),pos = list("b", "b"),
            height = unit(0.7, "line"),size = "auto", abbrv = list(TRUE, FALSE)) +
  guides(color = guide_legend(ncol = 6)) + # Arrange legend items in 1 row (adjust 'ncol' as needed)
  proxy_theme

# Create a plot grid for CO2 (with its legend below), D13C, and D34S
figS4 <- plot_grid(
  CO2,   # CO2 plot (with legend below)
  D13C,  # D13C plot
  D34S,  # D34S plot
  ncol = 1,  # One column (stacked vertically)
  rel_heights = c(0.4, 0.3, 0.3),  # Make the CO2 plot shorter (0.4), the others are equal height (0.3)
  align = "v",  # Vertically align the plots
  axis = "lr"   # Ensure consistent axis across all plots (left-right)
)

# Save figure plots:
CairoPDF(file = file.path(save_dir, "Figure_S4.pdf"), width = 5, height = 7)
print(figS4)
dev.off()
ggsave(figS4,path = save_dir,filename = "Figure_S4.png", width = 5, height = 7, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure S5: Alaska Comparison plots ----
# -------------------------------------------------------------------------- #

setwd(main_dir)
ox = c(20.95)
varnames = c("pftalbiomass","burnedf") 
gplot=c()
ex <- c("both") 
fname <- "_both_oct.nc"
pname <- "oxygen_fire_productivity"

# loop over vars
for (v in 1:length(varnames)) {
  var <- varnames[v]
  
  tempox <- array(0, length(ox)) # set temporary array to hold global totals if needed
  
  oxfile <- (paste(LPJ_data,pname,"/",ox,fname,sep="")) # get correct file
  
  mdf <- decadalavg_mdf(oxfile,var)
  
  gplot[[v]] <- globalplot(mdf,diff=FALSE,sim = paste(ex), o2=ox, var=var) + O2_line_theme # create plot
}

alask_sim = gplot[[1]]
alask_ba = gplot[[2]]

brk=c(-1,0.01,0.2,0.5,1,1.5,2,2.5,3,6,9,12,15,18,21,24,27,Inf)
pal=c("white","#E6E9E2","#B2B6B2","#6E4A00","#A16706","#CB8200","#E0AB04","#FAD103","#FFF0A4",# colour scale for aboveground biomass:
      "#D6ED00","#A2DE03","#7FC108","#569A06","#308100","#006400","#014909","#080A04")
unit = "Aboveground \nbiomass \n(kgm-2)"
leglab = c("0-0.01","0.01-0.2","0.2-0.5","0.5-1","1-1.5","1.5-2","2-2.5","2.5-3","3-6","6-9","9-12",
           "12-15","15-18","18-21","21-24", "24-27", "27+")

# load observations
obs = Alaska_data

# Apply the breaks to categorize the AGB values
obs$breaks <- cut(obs$Average_AGB/10, breaks = brk, include.lowest = TRUE, right = FALSE)
obs$breaks <- factor(obs$breaks)  # automatically drops unused levels

alask_obs = ggplot() +
  geom_tile(data = obs, aes(x = Lon, y = Lat, fill = breaks)) +
  scale_fill_manual(values = pal, drop = FALSE, labels = leglab) +  # Use custom colors
  xlim(range(obs$Lon)) + ylim(range(obs$Lat)) +
  guides(fill = guide_legend(direction = "vertical", title.position = "top",title = paste(unit)))+
  coord_equal(ylim = c())+
  labs(title = paste("", sep=""),
       tag = paste("B"),
       x= "Longitude",
       y= "Latitude")+
  alaska_theme

# Extract the legend from one of the plots (assuming both plots share the same legend)
legend <- get_legend(alask_obs)

lat_range <- c(55, 72)
lon_range <- c(-170, -130)

# Set consistent limits
alask_sim <- alask_sim + coord_equal(xlim = lon_range, ylim = lat_range) +
  labs(title = paste("", sep=""),
       tag = paste("A"),
       x= "Longitude",
       y= "Latitude") +
  guides(fill = guide_legend(direction = "vertical", title.position = "top"))+
  alaska_theme

alask_ba <- alask_ba + coord_equal(xlim = lon_range, ylim = lat_range) +
  labs(title = paste("", sep=""),
       tag = paste("C"),
       x= "Longitude",
       y= "Latitude") +
  guides(fill = guide_legend(direction = "vertical", title.position = "top"))+
  alaska_theme

alask_obs <- alask_obs + coord_equal(xlim = lon_range, ylim = lat_range)
legend2 <- get_legend(alask_ba)

figS5 <- grid.arrange(
  # Left column: three plots stacked
  arrangeGrob(
    alask_sim + theme(legend.position = "none"), 
    alask_obs + theme(legend.position = "none"), 
    alask_ba + theme(legend.position = "none"), 
    nrow = 3),                 # Stack three plots vertically
  
  # Right column: first legend spanning two rows for the first two plots, second legend for the third plot
  arrangeGrob(
    legend,                    # Corresponds to the first two plots, spanning two rows
    nullGrob(),                # Empty row to maintain spacing
    legend2,                   # Corresponds to the third plot
    heights = c(2, 0.5, 1),    # Adjust heights to balance legends and empty row
    nrow = 3),                 # Three rows to match the plot rows
  
  # Define layout: two columns (plots & legends)
  ncol = 2, widths = c(10, 2)  # Adjust the widths of the columns
)

# Save figure plots:
ggsave(figS5,path = save_dir,filename = "Figure_S5.pdf", width = 5.5, height = 6.5, units = c("in"))
ggsave(figS5,path = save_dir,filename = "Figure_S5.png", width = 5.5, height = 6.5, units = c("in"),dpi=900)

# -------------------------------------------------------------------------- #
## Figure S6: Global Tree cover under high CO2 ----
# -------------------------------------------------------------------------- #

figS6 <- plot_oxygen_grid(
  experiments = c("fire","product","both"),
  folders     = c("fire_only","photo_only","fire_and_photo"),
  varname     = "treecov",
  oxygen_levels = c(20.95, 25, 35),
  LPJ_data    = paste0(LPJ_data,"late_cret_2025/"),
  code_dir    = code_dir,
  climate_config = "LPJ_high_CO2"
)

# Save figure plots:
ggsave(figS6,path = save_dir,filename = "Figure_S6.pdf", width = 7.24, height = 4.6, units = c("in"))
ggsave(figS6,path = save_dir,filename = "Figure_S6.png", width = 7.24, height = 4.6, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure S7: Global Tree cover under high Temp  ----
# -------------------------------------------------------------------------- #

figS7 <- plot_oxygen_grid(
  experiments = c("fire","product","both"),
  folders     = c("fire_only","photo_only","fire_and_photo"),
  varname     = "treecov",
  oxygen_levels = c(20.95, 25, 35),
  LPJ_data    = paste0(LPJ_data,"late_cret_2025/"),
  code_dir    = code_dir,
  climate_config = "LPJ_high_temp"
)

# Save figure plots:
ggsave(figS7,path = save_dir,filename = "Figure_S7.pdf", width = 7.24, height = 4.6, units = c("in"))
ggsave(figS7,path = save_dir,filename = "Figure_S7.png", width = 7.24, height = 4.6, units = c("in"),dpi=900)

# -------------------------------------------------------------------------- #
## Figure S8: Global Tree cover under high CO2 & Temp  ----
# -------------------------------------------------------------------------- #

figS8 <- plot_oxygen_grid(
  experiments = c("fire","product","both"),
  folders     = c("fire_only","photo_only","fire_and_photo"),
  varname     = "treecov",
  oxygen_levels = c(20.95, 25, 35),
  LPJ_data    = paste0(LPJ_data,"late_cret_2025/"),
  code_dir    = code_dir,
  climate_config = "LPJ_high_CO2_temp"
)

# Save figure plots:
ggsave(figS8,path = save_dir,filename = "Figure_S8.pdf", width = 7.24, height = 4.6, units = c("in"))
ggsave(figS8,path = save_dir,filename = "Figure_S8.png", width = 7.24, height = 4.6, units = c("in"),dpi=900)


# -------------------------------------------------------------------------- #
## Figure S9: Global total vegetation over O2 under different climates ----
# -------------------------------------------------------------------------- #

# List of climate configs and short names for supplementary figure
supp_climate_configs <- c("LPJ_high_CO2", "LPJ_high_temp", 
                          "LPJ_high_CO2_temp")
short_names <- c("high_CO2", "high_temp", "high_CO2_temp")

# Here: biomass (A) and forest cover (B), (a)-(c) for climate configs
supp_labels <- c("A(a)", "B(a)", "A(b)", 
                 "B(b)", "A(c)", "B(c)")

# Store all plots and legends
supp_vplots <- list()
supp_legends <- list()

for (i in seq_along(supp_climate_configs)) {
  tmp <- make_climate_oxline_plots(
    climate_config = supp_climate_configs[i],
    short_name = short_names[i],
    O2_line_theme = O2_line_theme,
    base_dir = data_dir,
    custom_labels = supp_labels[((i-1)*2 + 1):((i-1)*2 + 2)] # pass 2 labels per climate
  )
  
  supp_vplots[[i]] <- tmp$vplot
  supp_legends[[i]] <- tmp$legend
}

# Extract a single legend (assuming all plots share the same style)
legend_supp <- supp_legends[[1]]

# Build the grid 
supp_clim_oxline <- cowplot::plot_grid(
  supp_vplots[[1]][[1]] + theme(legend.position = "none"),
  supp_vplots[[1]][[2]] + theme(legend.position = "none"),
  supp_vplots[[2]][[1]] + theme(legend.position = "none"),
  supp_vplots[[2]][[2]] + theme(legend.position = "none"),
  supp_vplots[[3]][[1]] + theme(legend.position = "none"),
  supp_vplots[[3]][[2]] + theme(legend.position = "none"),
  nrow = 3, align = "hv"
)

# Add legend underneath
fig_s9 <- cowplot::plot_grid(supp_clim_oxline, legend_supp, ncol = 1, rel_heights = c(1, 0.08))

print(fig_s9)

# Save supplementary figure
ggsave(fig_s9, path = save_dir, filename = "Figure_S9.pdf", width = 5.2, height = 6.5, units = "in")
ggsave(fig_s9, path = save_dir, filename = "Figure_S9.png", width = 5.2, height = 6.5, units = "in", dpi = 900)


# -------------------------------------------------------------------------- #
## Figure S10: Total vegetation over O2 under different climates ----
# -------------------------------------------------------------------------- #

# Setup variables
climate_configs <- c("LPJ_high_CO2", "LPJ_high_temp", "LPJ_high_CO2_temp","LPJ_high_CO2_temp_precip")
short_names <- c("high_CO2", "high_temp", "high_CO2_temp", "high_CO2_temp_precip")

experiments <- c("FIRE_AND_PHOTO","FIRE_ONLY", "PHOTO_ONLY")
exp_labels <- c( "combined","fire", "photorespiration")
line_colors <- c("#3A4454", "#F15025", "#1A8FE3")

# Normalized biomass column suffix
vval <- "pftalbiomass_norm"
ylab_text <- "terrestrial biomass (normalised)"
ylim_val <- c(0, 4)

# Prepare a list to hold the plots
plots <- list()

# Prediction oxygen values
preddat <- data.frame(Ox = seq(17, 35, by = 0.5))  # oxygen range

for (cc in seq_along(climate_configs)) {
  config <- climate_configs[cc]
  short <- short_names[cc]
  
  # Load data
  sheet_path <- paste0(data_dir,"/LPJLMfire_output/totals/climate_configurations/", 
                       config, "_master.xlsx")
  my_data <- read_excel(sheet_path)
  my_data <- as.data.frame(my_data)
  
  # Prepare matrix to hold fits
  splines <- matrix(ncol = 4, nrow = nrow(preddat))
  
  # Fit GAM for each experiment
  for (exp in seq_along(experiments)) {
    colname <- paste0(experiments[exp], "_", vval)
    exdf <- data.frame(Ox = my_data$O2, var = as.numeric(my_data[[colname]]))
    
    bsmod <- gam(var ~ s(Ox, k = 14, bs = "bs", m = c(3, 2)), data = exdf, method = "REML")
    pred <- predict(bsmod, preddat)
    
    splines[, exp] <- pred
  }
  
  # Prepare long format for plotting
  plotdat <- data.frame(
    Ox = rep(preddat$Ox, times = 3),
    fit = c(splines[,1], splines[,2], splines[,3]),
    Experiment = factor(rep(exp_labels, each = nrow(preddat)), levels = exp_labels)
  )
  
  # Plot only smooth lines, no ribbons
  p <- ggplot(plotdat, aes(x = Ox, y = fit, color = Experiment)) +
    geom_line(size = 0.6) +
    scale_color_manual(values = line_colors, name = "") +
    scale_x_continuous(limits = c(17, 35), breaks = seq(17, 35, 2)) +
    # scale_y_continuous(limits = ylim_val) +
    expand_limits(y = 0) + 
    labs(x = "Atmospheric oxygen concentration (%)", y = ylab_text,
         # title = paste0(climate_labels[cc])) +
         title = "") +
    O2_line_theme
  
  plots[[cc]] <- p
}

# Extract legend from first plot
leg <- g_legend(plots[[1]] + theme(legend.position = "bottom", legend.direction = "horizontal"))


# Remove legends from individual plots
plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))

# Stack vertically with labels and add legend at bottom
fig_s10 <- plot_grid(
  plot_grid(plotlist = plots_no_legend, ncol = 1, labels = c("A", "B", "C","D"), align = "v"),
  leg,
  ncol = 1,
  rel_heights = c(1, 0.025)
)

print(fig_s10)

# Save supplementary figure
ggsave(fig_s10, path = save_dir, filename = "Figure_S10.pdf", width = 5.2, height = 6.5, units = "in")
ggsave(fig_s10, path = save_dir, filename = "Figure_S10.png", width = 5.2, height = 6.5, units = "in", dpi = 900)

