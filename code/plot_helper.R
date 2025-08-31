## Script of helper functions and settings for plots

# --- Function to save a legend seperately from a plot --- #
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# --- Function to load and process data --- #
load_and_melt_data <- function(file_path, cols) {
  data <- read_excel(file_path)
  melted_data <- melt(data.frame(data[, cols]), id.vars = "O2_conc")
  return(melted_data)
}

# --- Function to compute min and max values and fit loess models --- #
compute_min_max_and_loess <- function(melted_data, span = 0.5) {
  min_max_data <- melted_data %>%
    group_by(O2_conc) %>%
    summarise(
      min_value = min(value),
      max_value = max(value)
    )
  min_loess <- loess(min_value ~ O2_conc, data = min_max_data, span = span)
  max_loess <- loess(max_value ~ O2_conc, data = min_max_data, span = span)
  return(list(min_max_data = min_max_data, min_loess = min_loess, max_loess = max_loess))
}

# --- Function to generate smoothed values and create a data frame --- #
generate_smoothed_data <- function(loess_models, x_seq) {
  min_smooth <- predict(loess_models$min_loess, x_seq)
  max_smooth <- predict(loess_models$max_loess, x_seq)
  smoothed_data <- data.frame(O2_conc = x_seq, min_smooth = min_smooth, max_smooth = max_smooth)
  return(smoothed_data)
}

# --- theme definition for O2 line plots ---#
O2_line_theme <- 
  theme(  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.margin = margin(0.1,0.1,0.05,0.1,unit = "cm"),
    axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0)),
    axis.title = element_text(size = 8,vjust=-5),
    panel.border = element_rect(fill=NA, colour="black", linewidth=0.1),
    plot.title = element_text(size=9,face="bold"),
  )

# --- Theme definition for COPSE plots --- #
copse_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8, colour = "grey30"),
        axis.title = element_text(size = 8, vjust=-5),
        plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"),
        plot.tag = element_text(face = "bold", colour = "black"),
        plot.tag.position = "topleft",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.1),
        legend.key.height = unit(1, "lines"),
        legend.key.size =  unit(0.5, "cm"),
        legend.text  = element_text(size = 7),
        legend.box.spacing = unit(0.5, "cm"),
        legend.background = element_rect(color = "grey80", linewidth = 0.2)
  )

# --- Theme for Alaska plots --- #
alaska_theme <- theme(  
  legend.key.size = unit(0.4, "cm"),
  legend.key = element_rect(colour = '#bdbdbd', size = 0.1),
  legend.title = element_text(size = 7),
  legend.position = "right",
  legend.box.margin = margin(0,0,0,0, unit="cm"),
  panel.background = element_rect(fill = "white"),
  legend.text = element_text(size=7),
  legend.box.background = element_rect(colour = "grey33", size = 0.1),
  legend.margin = margin(0,0,0,0,unit = "cm"),
  panel.border = element_rect(fill=NA, colour="black", size=0.1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_blank(),
  axis.text =  element_text(size=8),
  axis.ticks = element_blank(),
  axis.title = element_text(size=8),
  plot.title = element_text(size=9,hjust = 0, vjust = 0,margin = margin(0,0,2,0)),
  plot.tag = element_text(size = 9, face="bold",hjust = -0.5, vjust = -1.),
  plot.tag.position = c(0, 1),
  plot.margin = unit(c(0.5,0,0,0), "cm")
)

# --- Theme for proxy plot ---
proxy_theme <- theme_classic()+
  theme(legend.position = "bottom",            # Legend at the bottom of CO2 plot
        legend.key.size = unit(0.5, "cm"),     # Make legend key smaller
        legend.text = element_text(size = 7),  # Smaller legend text
        legend.title = element_text(size = 7), # Smaller legend title
        axis.text = element_text(size = 8),    # Smaller axis text
        axis.title = element_text(size = 8))   # Smaller axis title

# --- Color and linetype mappings for COPSE plots --- #
copse_colors <- c("Fire & \nphotorespiration"="black",
                  "No vegetation-based \noxygen feedback"="gray50",
                  "Photorespiration only \n(no fire)"="#1A8FE3",
                  "Fire only \n(no photorespiration)"="#F15025",
                  "Combined"="#3A4454")

copse_linetypes <- c("Fire & \nphotorespiration" = "dashed",
                     "Fire only \n(no photorespiration)"="solid",
                     "Photorespiration only \n(no fire)"="solid",
                     "Combined"="solid",
                     "No vegetation-based \noxygen feedback" = "dotted",
                     "Mills"="dotdash")

# --- Function to create a 3x3 oxygen × scenario grid plot for LPJ-LMfire data --- #
plot_oxygen_grid <- function(experiments,
                             folders,
                             varname,
                             oxygen_levels,
                             LPJ_data,
                             code_dir,
                             climate_config = NULL,   # <-- only one at a time
                             lab    = c("A","B","C",""),
                             lablab = c("(a)","(b)","(c)")) {
  
  gplot <- list()
  gp <- 1
  
  for (exp in seq_along(experiments)) {
    ex <- experiments[exp]
    folder <- folders[exp]
    
    # --- Case 1: Figure 1 (no climate_config) ---
    if (is.null(climate_config)) {
      if (ex == "fire") {
        fname <- "_fire_july22.nc"
        sim_label <- "fire"
      } else if (ex == "product") {
        fname <- "_andre.nc"
        sim_label <- "photorespiration"
      } else if (ex == "both") {
        fname <- "_both_oct.nc"
        sim_label <- "combined"
      }
    }
    
    # --- Case 2: Climate config (Figure 3 etc.) ---
    if (!is.null(climate_config)) {
      if (ex == "fire") {
        pname <- "fire_only"; sim_label <- "fire"
      } else if (ex == "product") {
        pname <- "photo_only"; sim_label <- "photorespiration"
      } else if (ex == "both") {
        pname <- "fire_and_photo"; sim_label <- "combined"
      }
    }
    
    for (o in seq_along(oxygen_levels)) {
      ox <- oxygen_levels[o]
      
      # Build nc file path
      if (is.null(climate_config)) {
        nc_file <- file.path(LPJ_data, folder, paste0(ox, fname))
      } else {
        fname <- paste0(climate_config, "_", pname, "_O2", ox, ".nc")
        nc_file <- file.path(LPJ_data, climate_config, pname, fname)
      }
      
      if (!file.exists(nc_file)) next
      
      # Special case: tree cover stored as forestcov
      if (varname == "treecov") {
        mdf <- decadalavg_mdf(nc_file, var = "forestcov", path = code_dir)
      } else {
        mdf <- decadalavg_mdf(nc_file, var = varname, path = code_dir)
      }
      
      # Labels (A(i), B(ii), etc.)
      labs <- paste0(lab[exp], lablab[o])
      
      gplot[[gp]] <- globalplot(
        mdf,
        diff = FALSE,
        sim = sim_label,
        plotlab = labs,
        o2 = ox,
        var = varname
      )
      gp <- gp + 1
    }
  }
  
  if (length(gplot) < 2) {
    stop("No plots were generated. Check that the .nc files exist and paths are correct.")
  }
  
  # Grab legend from second plot
  leg <- g_legend(gplot[[2]])
  
  # Apply common theme
  gplots_no_legend <- lapply(gplot, function(p) {
    p + theme(
      legend.position = "none",
      plot.margin = unit(c(-0.5, 0.1, -1, -0.1), "cm")
    )
  })
  
  # Arrange in 3x3 grid
  grid_plot <- grid.arrange(
    arrangeGrob(
      grobs = gplots_no_legend,
      nrow = 3
    ),
    leg,
    heights = c(10, 1.3),
    top = ""
  )
  
  return(grid_plot)
}


# --- Function to create line plots over oxygen for a given climate 
make_climate_oxline_plots <- function(climate_config, short_name, 
                                      experiments = c("FIRE_ONLY", "PHOTO_ONLY", "FIRE_AND_PHOTO"),
                                      exp_labels = c("fire", "photorespiration", "combined"), 
                                      line_colors = c("#3A4454", "#F15025", "#1A8FE3"),
                                      vval = c("pftalbiomass_val", "forestcov_val"),
                                      ylabs = c("aboveground biomass (PgC)", expression("forest cover (km"^2*" x10"^6*")")),
                                      ylims = c(1750, 100),
                                      var_letters = c("A", "B"),
                                      lab_labels = c("(a)", "(b)", "(c)"),
                                      O2_line_theme,
                                      base_dir,
                                      custom_labels = NULL) {   # optional override
  
  vplot <- list()
  
  # Loop over variables: 1 = biomass, 2 = forest cover
  for (var in 1:2) {
    sheet_path <- paste0(base_dir, "/LPJLMfire_output/totals/climate_configurations/", climate_config, "_master.xlsx")
    
    my_data <- readxl::read_excel(sheet_path)
    my_data <- as.data.frame(my_data)
    
    preddat <- data.frame(Ox = seq(21, 35, by = 0.5))
    splines <- matrix(ncol = length(experiments), nrow = nrow(preddat))
    
    for (exp in seq_along(experiments)) {
      colname <- paste0(experiments[exp], "_", vval[var])
      exdf <- data.frame(Ox = my_data$O2, var = as.numeric(my_data[[colname]]))
      
      bsmod <- mgcv::gam(var ~ s(Ox, k=14, bs="bs", m=c(3,2)), data=exdf, method="REML")
      pred <- predict(bsmod, preddat, se.fit = TRUE)
      
      if (var == 2) splines[, exp] <- pred$fit / 1e6 else splines[, exp] <- pred$fit
    }
    
    # Determine plot title
    plot_title <- if (!is.null(custom_labels)) custom_labels[var] else var_letters[var]
    
    # Build ggplot
    p <- ggplot() +
      geom_hline(yintercept = 0, col = "grey60", linetype = "dotted")
    
    for (exp in seq_along(experiments)) {
      p <- p + geom_line(aes_(x = preddat$Ox, y = splines[, exp],
                              color = exp_labels[exp],
                              linetype = exp_labels[exp]),
                         size = 0.4)
    }
    
    p <- p +
      scale_color_manual(values = line_colors, name = "") +
      scale_linetype_manual(values = rep("solid", length(experiments)), name = "") +
      scale_x_continuous(limits = c(20.95, 35), breaks = seq(20, 36, 2)) +
      scale_y_continuous(limits = c(0, ylims[var]), breaks = seq(0, ylims[var], ylims[var]/10)) +
      labs(x = "Oxygen (% of atmosphere)", y = ylabs[[var]], title = plot_title) +
      O2_line_theme
    
    vplot[[var]] <- p
  }
  
  # Extract legend
  leg <- g_legend(vplot[[1]])
  
  return(list(vplot = vplot, legend = leg))
}
