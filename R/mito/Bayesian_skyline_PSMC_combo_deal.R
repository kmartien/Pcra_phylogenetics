# ==============================================================================
# COMBINED DEMOGRAPHIC ANALYSIS: PSMC/hPSMC (Nuclear) & Bayesian Skyline (mtDNA)
# ==============================================================================
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
load('data/broad_stratum_colors.rda')

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION & PATHS
# ------------------------------------------------------------------------------
PSMC_FOLDER     <- "hPSMC_Dsuite/hPSMC"
BSP_FILE_PATH   <- "BEAST/xml/aln3-all/aln3-all_skyline_reconstruction" 
# Biological Parameters (Hernandez et al. 2026)
GENERATION_TIME <- 24 
SEX_RATIO_SCALE <- 2.0 
TIME_SCALE      <- 1e6 

# Shared Axis Adjustments
XMIN <- 0
XMAX <- 225000

# Primary Y-Axis Limits (Nuclear PSMC Linear Scale - Now True Absolute Ne)
YMIN_PSMC <- 0
YMAX_PSMC <- 180000  # Shifted from 18 to 180,000 for true scale plotting

# ------------------------------------------------------------------------------
# 2. LOAD AND PROCESS BAYESIAN SKYLINE DATA (MITOCHONDRIAL)
# ------------------------------------------------------------------------------
skyline_raw <- read_tsv(BSP_FILE_PATH, skip = 1)

skyline_data <- skyline_raw %>%
  mutate(
    Ne_median = (median * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME,
    Ne_upper  = (upper * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME,
    Ne_lower  = (lower * SEX_RATIO_SCALE * TIME_SCALE) / GENERATION_TIME,
    time_scaled = time * TIME_SCALE
  ) %>%
  filter(time_scaled >= XMIN & time_scaled <= XMAX)

# Automatically capture the data boundaries in log10 space to scale perfectly
BSP_LOG_MIN <- floor(log10(min(skyline_data$Ne_lower)))
BSP_LOG_MAX <- ceiling(log10(max(skyline_data$Ne_upper)))

# Transform log-scale BSP data into the linear 0-180,000 PSMC axis space
skyline_data <- skyline_data %>%
  mutate(
    Ne_median_scaled = (log10(Ne_median) - BSP_LOG_MIN) / (BSP_LOG_MAX - BSP_LOG_MIN) * YMAX_PSMC,
    Ne_upper_scaled  = (log10(Ne_upper) - BSP_LOG_MIN) / (BSP_LOG_MAX - BSP_LOG_MIN) * YMAX_PSMC,
    Ne_lower_scaled  = (log10(Ne_lower) - BSP_LOG_MIN) / (BSP_LOG_MAX - BSP_LOG_MIN) * YMAX_PSMC
  )

# ------------------------------------------------------------------------------
# 3. LOAD AND PROCESS PSMC / hPSMC DATA (NUCLEAR)
# ------------------------------------------------------------------------------
allfiles <- list.files(path = PSMC_FOLDER, pattern = "txt")

pop_metadata <- data.frame(
  id         = c("Pcra_z0018462", "Pcra_z0027510", "Pcra_z0045928", 
                 "hPSMC_z0018462_z0027510", "hPSMC_z0018462_z0045928", "hPSMC_z0027510_z0045928"),
  clean_name = c("ETP", "ATL", "MHI", "hPSMC:ETP-ATL", "hPSMC:ETP-MHI", "hPSMC:ATL-MHI"),
  stringsAsFactors = FALSE
)

psmc_main_list <- list()
psmc_boot_list <- list()

for (i in 1:nrow(pop_metadata)) {
  pid   <- pop_metadata$id[i]
  pname <- pop_metadata$clean_name[i]
  
  p_files <- allfiles[grep(pattern = pid, allfiles)]
  if (length(p_files) == 0) next
  
  # Main trajectories (Scaled up to absolute Ne values)
  main_df <- read.table(file.path(PSMC_FOLDER, p_files[1])) %>%
    rename(time = V1, Ne = V2) %>%
    mutate(Ne = Ne * 10000) %>% 
    filter(time >= XMIN & time <= XMAX) %>%
    mutate(pop = pname)
  psmc_main_list[[pid]] = main_df
  
  # Bootstrap replicates (Scaled up to absolute Ne values; skip for hPSMC files)
  if (length(p_files) > 1 & !grepl("hPSMC", pid)) {
    for (j in 2:length(p_files)) {
      boot_df <- read.table(file.path(PSMC_FOLDER, p_files[j])) %>%
        rename(time = V1, Ne = V2) %>%
        mutate(Ne = Ne * 10000) %>% 
        filter(time >= XMIN & time <= XMAX) %>%
        mutate(pop = pname, rep = j)
      psmc_boot_list[[paste0(pid, "_", j)]] <- boot_df
    }
  }
}

psmc_main_data <- bind_rows(psmc_main_list)
psmc_boot_data <- bind_rows(psmc_boot_list)

# ------------------------------------------------------------------------------
# 4. COLOR, LINETYPE, AND UNIFIED LEGEND CONFIGURATION
# ------------------------------------------------------------------------------
LEGEND_NAME   <- "Lineage / Analysis"
LEGEND_BREAKS <- c("ETP", "ATL", "MHI", "hPSMC:ETP-ATL", "hPSMC:ETP-MHI", "hPSMC:ATL-MHI", "Mitochondrial BSP (Median)")

plot_colors <- c(
  "ETP"                        = rgb(229, 198, 27, maxColorValue = 255),
  "ATL"                        = rgb(72, 21, 104, maxColorValue = 255),
  "MHI"                        = rgb(117, 117, 117, maxColorValue = 255),
  "hPSMC:ETP-ATL"              = "black",
  "hPSMC:ETP-MHI"              = "black",
  "hPSMC:ATL-MHI"              = "black",
  "Mitochondrial BSP (Median)" = "blue"
)

plot_linetypes <- c(
  "ETP"                        = "solid",
  "ATL"                        = "solid",
  "MHI"                        = "solid",
  "hPSMC:ETP-ATL"              = "solid",
  "hPSMC:ETP-MHI"              = "dashed",
  "hPSMC:ATL-MHI"              = "dotted",
  "Mitochondrial BSP (Median)" = "solid"
)

# Set up breaks and clean labels for the log secondary axis
bsp_breaks <- BSP_LOG_MIN:BSP_LOG_MAX
bsp_labels <- prettyNum(10^bsp_breaks, big.mark = ",")

# ------------------------------------------------------------------------------
# 5. GENERATE DUAL-AXIS GGPLOT
# ------------------------------------------------------------------------------
p <- ggplot() +
  # Place mitochondrial HPD ribbon in background first
  geom_ribbon(data = skyline_data, 
              aes(x = time_scaled, ymin = Ne_lower_scaled, ymax = Ne_upper_scaled), 
              fill = "skyblue", alpha = 0.25) +
  
  # Mitochondrial Median Line (Maps both color & linetype inside aes to unify legends)
  geom_line(data = skyline_data, 
            aes(x = time_scaled, y = Ne_median_scaled, 
                color = "Mitochondrial BSP (Median)", 
                linetype = "Mitochondrial BSP (Median)"), 
            linewidth = 1) +
  
  # Nuclear Bootstrap Replicates (Thin & Translucent)
  geom_step(data = psmc_boot_data, 
            aes(x = time, y = Ne, group = interaction(pop, rep), color = pop), 
            alpha = 0.04, linewidth = 0.3, direction = "hv", show.legend = FALSE) +
  
  # Nuclear Main Trajectories & hPSMC
  geom_step(data = psmc_main_data, 
            aes(x = time, y = Ne, color = pop, linetype = pop), 
            linewidth = 0.9, direction = "hv") +
  
  # Formatting Axes
  scale_x_continuous(
    limits = c(XMIN, XMAX),
    breaks = seq(XMIN, XMAX, by = 50000),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(YMIN_PSMC, YMAX_PSMC),
    breaks = seq(YMIN_PSMC, YMAX_PSMC, by = 20000), # Ticks every 20,000 individuals
    labels = scales::label_comma(),
    name = expression("Nuclear Effective Population Size (" * italic(N)[e] * ") [Linear Scale]"),
    
    # Linear transformation vector that maps grid-lines back to true Log-values
    sec.axis = sec_axis(
      trans = ~ BSP_LOG_MIN + . * (BSP_LOG_MAX - BSP_LOG_MIN) / YMAX_PSMC,
      name = expression("Mitogenome Effective Population Size (" * italic(N)[e] * ") [Log Scale]"),
      breaks = bsp_breaks,
      labels = bsp_labels
    )
  ) +
  
  # Merging scale setups using identical Names and Breaks layouts
  scale_color_manual(name = LEGEND_NAME, breaks = LEGEND_BREAKS, values = plot_colors) +
  scale_linetype_manual(name = LEGEND_NAME, breaks = LEGEND_BREAKS, values = plot_linetypes) +
  
  labs(x = "Years before present") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    axis.title.y.right = element_text(color = "blue", margin = margin(l = 10)),
    axis.title.y.left  = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")
  )

# ------------------------------------------------------------------------------
# 6. EXPORT GRAPHIC
# ------------------------------------------------------------------------------
ggsave("results-raw/combined_nuclear_mitogenome_demography.pdf", 
       plot = p, width = 8.5, height = 5, dpi = 300)