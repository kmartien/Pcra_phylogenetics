
# Install and load required packages if you haven't already
# install.packages("ggplot2")
# install.packages("readr")
# install.packages("scales")
library(ggplot2)
library(readr)
library(scales) 

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
FILE_PATH       <- "BEAST/xml/aln3-all/aln3-all_skyline_reconstruction" 
GENERATION_TIME <- 24.0 #from Hernandez et al. 2026                                
SEX_RATIO_SCALE <- 2.0                                
TIME_SCALE      <- 1000000                            # Scaled to 1e6 to convert decimals to years

# ==============================================================================
# 2. LOAD AND CONVERT DATA
# ==============================================================================
skyline_data <- read_tsv(FILE_PATH, skip = 1)

# Apply demographic conversion formulas
skyline_data$Ne_median <- (skyline_data$median * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_upper  <- (skyline_data$upper * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_lower  <- (skyline_data$lower * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME

# Scale time to single years ago
skyline_data$time_scaled <- skyline_data$time * TIME_SCALE

# CRITICAL FOR LOG-SCALE: Remove the 0.0 timepoint to prevent log(0) errors
skyline_data <- skyline_data[skyline_data$time_scaled > 0, ]

# ==============================================================================
# 3. GENERATE PLOT
# ==============================================================================
p <- ggplot(skyline_data, aes(x = time_scaled)) +
  # Draw the 95% HPD interval as a shaded ribbon
  geom_ribbon(aes(ymin = Ne_lower, ymax = Ne_upper), 
              fill = "skyblue", alpha = 0.4) +
  
  # Draw the median population trajectory line
  geom_line(aes(y = Ne_median), 
            color = "blue", linewidth = 1) +
  
  # Linear scale for Y-axis (Standard numbers, no log compression)
  scale_y_continuous(labels = comma) +
  
  # Log scale for X-axis (To match PSMC timeline spacing)
  scale_x_log10(labels = comma, 
                breaks = 10^(1:7)) + # Major breaks at 10, 100, 1k, 10k, 100k, 1M...
  
  # Adds log ticks strictly to the bottom (X-axis)
  annotation_logticks(sides = "b", color = "gray60") + 
  
  # Labels and clean styling
  labs(
    x = "Time (Years Ago)",
    y = expression("Total Effective Population Size (" * italic(N)[e] * ")"),
    title = "Historical Demographic Trajectory",
    subtitle = "Coalescent Bayesian Skyline Plot (Linear-Log Scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5)
  )

# Save the plot directly to high-res PDF or PNG
ggsave("R/mito/bayesian_skyline_linearY_logX_plot.pdf", width = 8, height = 5, dpi = 300)
