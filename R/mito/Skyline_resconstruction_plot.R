# Install and load required packages if you haven't already
# install.packages("ggplot2")
# install.packages("readr")
library(ggplot2)
library(readr)

# ==============================================================================
# 1. USER CONFIGURATION (Adjust these based on your organism)
# ==============================================================================
# Relative path from the main project folder root to your data file
FILE_PATH       <- "BEAST/xml/aln3-all/aln3-all_skyline_reconstruction" 
GENERATION_TIME <- 24 # Generation time (in years); from Hernandez et al. 2026
SEX_RATIO_SCALE <- 2.0 # 2.0 scales mtDNA to total Ne (males+females)
TIME_SCALE      <- 1e6 # Change to 1e6 if X-axis is in millions of years

# ==============================================================================
# 2. LOAD AND CONVERT DATA
# ==============================================================================
# We use skip = 1 to bypass the metadata string Tracer inserts at line 1
skyline_data <- read_tsv(FILE_PATH, skip = 1)

# Apply the mathematical demographic conversion directly to the dataframe
skyline_data$Ne_median <- (skyline_data$median * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_upper  <- (skyline_data$upper * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_lower  <- (skyline_data$lower * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME

# Scale time if necessary (e.g., converting to calendar years ago)
skyline_data$time_scaled <- skyline_data$time * TIME_SCALE

# ==============================================================================
# 3. GENERATE PUBLICATION-READY PLOT
# ==============================================================================
ggplot(skyline_data, aes(x = time_scaled)) +
  # Draw the 95% HPD interval as a shaded ribbon
  geom_ribbon(aes(ymin = Ne_lower, ymax = Ne_upper), 
              fill = "skyblue", alpha = 0.4) +
  
  # Draw the median population trajectory line
  geom_line(aes(y = Ne_median), 
            color = "blue", linewidth = 1) +
  
  # Log-transform the Y-axis
  scale_y_log10() +
  
  # Labels and clean styling
  labs(
    x = "Years before present",
    y = expression("Total Effective Population Size (" * italic(N)[e] * ")")
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 225000)) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5)
  )

# Optional: Save the plot directly to high-res PDF or PNG
# By default, ggsave saves to your working directory (project root). 
# You can change this path if you want the image saved elsewhere, e.g.:
 ggsave("results-raw/bayesian_skyline_plot.pdf", width = 8, height = 5, dpi = 300)