library(ggplot2)
library(readr)
library(scales)

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
FILE_PATH       <- "BEAST/xml/aln3-all/aln3-all_skyline_reconstruction" 
GENERATION_TIME <- 24 
SEX_RATIO_SCALE <- 2.0 
TIME_SCALE      <- 1e6 

# ==============================================================================
# 2. LOAD AND CONVERT DATA
# ==============================================================================
skyline_data <- read_tsv(FILE_PATH, skip = 1)

skyline_data$Ne_median <- (skyline_data$median * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_upper  <- (skyline_data$upper * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$Ne_lower  <- (skyline_data$lower * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME
skyline_data$time_scaled <- skyline_data$time * TIME_SCALE

# ==============================================================================
# 3. GENERATE IDENTICAL-THEME LOG PLOT
# ==============================================================================
ggplot(skyline_data, aes(x = time_scaled)) +
  geom_ribbon(aes(ymin = Ne_lower, ymax = Ne_upper), 
              fill = "skyblue", alpha = 0.4) +
  geom_line(aes(y = Ne_median), 
            color = "blue", linewidth = 1) +
  
  # Custom scale mapping: Forces explicit 50k blocks and regular plain integers
  scale_x_continuous(
    breaks = seq(0, 225000, by = 50000), 
    labels = scales::label_comma()
  ) +
  
  # Log scale configuration with explicit micro-breaks for tick mark generation
  scale_y_log10(
    breaks = c(
      seq(1000, 9000, by = 1000),
      seq(10000, 90000, by = 10000),
      seq(100000, 1000000, by = 100000)
    ),
    # Selectively labels major thresholds to prevent text crowding
    labels = function(x) {
      ifelse(x %in% c(1000, 5000, 10000, 50000, 100000, 500000, 1000000), scales::comma(x), "")
    }
  ) +
  
  labs(
    x = "Years before present",
    y = expression("Mitogenome Effective Population Size (" * italic(N)[e] * ")")
  ) +
  
  coord_cartesian(xlim = c(0, 225000)) +
  
  # Structural modifications to force matching design semantics
  theme_minimal(base_size = 11, base_family = "sans") + 
  theme(
    # Eliminate the outer frame lines entirely
    panel.border = element_blank(),
    
    # Establish standard matching original thin line (linewidth = 0.5)
    axis.line = element_line(color = "black", linewidth = 0.5),
    
    # Expose and style the y-axis ticks; keep x-axis tickless (linewidth = 0.5)
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.y = unit(4, "pt"),
    axis.ticks.x = element_blank(),
    
    # Grid modifications: Vertical only, styled thin light gray (linewidth = 0.5)
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")
  )

ggsave("results-raw/bayesian_skyline_plot.pdf", width = 6, height = 3.5, dpi = 300)