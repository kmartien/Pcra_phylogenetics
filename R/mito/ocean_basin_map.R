library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(maps)
library(mapdata)
library(sp)
library(scales)
library(patchwork)

source("R/mito/generate_ellipse_polygons.R")

# Custom function to format Longitude (handles the 360 wrap)
lon_formatter <- function(x) {
  # Convert values > 180 back to negative (e.g., 200 becomes -160)
  lon <- ifelse(x > 180, x - 360, x)
  # Assign E/W suffix
  suffix <- ifelse(lon < 0, "°W", ifelse(lon > 0, "°E", "°"))
  paste0(abs(lon), suffix)
}

# Custom function to format Latitude
lat_formatter <- function(x) {
  suffix <- ifelse(x < 0, "°S", ifelse(x > 0, "°N", "°"))
  paste0(abs(x), suffix)
}

load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/archive_data.rda')
load('../Pcra.database.data/data/ASPhyloCR.rda')
load('../Pcra.database.data/data/ASNGS.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

strata <- Pcra.strata |> 
  mutate(Broad = ifelse(Broad %in% c('Atl-weird', 'Atl-normal', 'Atl-unknown'), 'N_Atlantic', Broad),
  #       Broad = ifelse(Broad %in% c('MHI', 'NWHI'), 'HI_res', Broad),
         Broad = ifelse(Broad == 'Atl-South', 'S_Atlantic', Broad),
         Broad = ifelse(Broad == 'Spac', 'S_Pacific', Broad),
         Broad = ifelse(Broad == 'Wpac', 'NPac_Western', Broad),
         Broad = ifelse(Broad == 'Epac', 'NPac_Eastern', Broad),
         Broad = ifelse(Broad == 'CNP', 'NPac_Central', Broad),
         Broad = ifelse(Broad == 'Indian Ocean', 'Indian_Ocean', Broad))

samp.dat <- left_join(strata,
                      (select(archive_data, c(LABID,Latitude, Longitude)) |> 
                         rename(Animal.ID = LABID))) |> 
  left_join(CR.haps) |> 
  left_join(mito.haps) |> 
  mutate(
    analysis_set = ifelse(Animal.ID %in% ASNGS$Animal.ID, 'mito', 
                          ifelse(Animal.ID %in% ASPhyloCR$Animal.ID, 'CR', NA)),
    Longitude = ifelse(Longitude < 0,
                        ifelse(
                          Broad %in% c("NPac_Eastern", "NPac_Central", "S_Pacific", "MHI", "NWHI"), 
                          Longitude+360, 
                          Longitude),
                        Longitude)
      ) |> 
  filter(!is.na(analysis_set))

# Define desired stratum order
basin_order <- c(
  "N_Atlantic", 
  "S_Atlantic", 
  "Indian_Ocean", 
  "S_Pacific", 
  "NPac_Western", 
  "NPac_Central", 
  "NPac_Eastern",
  "Aus_resident",
  "MHI",
  "NWHI"
)

# Create a named vector for the "pretty" stratum labels
label_data <- tribble(
  ~Broad,          ~Mean_Long,  ~Mean_Lat,  ~Label_Name,
  "N_Atlantic",    -45,    40,  "North\nAtlantic",
  "S_Atlantic",    -20,   -25,  "South\nAtlantic",
  "Indian_Ocean",   60,     -7,  "Indian\nOcean",
  "S_Pacific",     230,   -10,  "South\nPacific",
  "NPac_Western",  96,    40,  "Western\nN. Pacific",
  "NPac_Central",  197.05908,    45,  "Central\nN. Pacific",
  "NPac_Eastern",  272,    41,  "Eastern\nN. Pacific",
  "Aus_resident",  86,   -29,  "Australia\nResident"#,
#  "HI_res",        202.30489,    20.74375,  "MHI & NWHI\nResidents"
)

# Create a named vector for the "pretty" stratum labels
placenames <- tribble(
  ~Name,          ~Mean_Long,  ~Mean_Lat,  ~Label_Name, ~Hjust,
  "Gabon",    20,    3,  "Gabon", 0,
  "Oman",    40,   25,  "Oman", 1
)

# Update samp.dat to make 'Broad' a factor with those levels
samp.dat <- samp.dat |>
  # this filters out two Korean market samples and a HI sample that could be
  # MHI resident or NWHI resident
  filter(Broad != 'Unknown') |> 
  # Move HI_res samples to the end of the data frame so they plot last
#  arrange(Broad == "HI_res") |> 
  arrange(analysis_set == "mito") |> 
  mutate(Broad = factor(Broad, levels = basin_order))

# Define initial ellipse parameters for each stratum
# center_x: Center Longitude (degrees)
#c enter_y: Center Latitude (degrees)
# a: Semi-major axis radius (ellipse width in degrees longitude)
# b: Semi-minor axis radius (ellipse height in degrees latitude)
# angle: Tilt rotation in degrees (counter-clockwise)
ellipse_params <- tribble(
  ~Broad,          ~center_x, ~center_y, ~a,   ~b,   ~angle,
  "N_Atlantic",    -67,       22,        30,   12,   -25,
  "S_Atlantic",    -23,      -28,        45,   12,    35,
  "Indian_Ocean",   92,       5,        45,   12,    -18,
  "S_Pacific",     215,      -15,        55,   12,   10,
  "NPac_Western",  141,       22,        30,   14,    -10,
  "NPac_Central",  202,       20,        27,   17,    0,
  "NPac_Eastern",  262,       19,        26,   10,    -33,
  "Aus_resident",  125,      -16,         8,    3,     24#,
#  "HI_res",        198,       21,         8,    4,    15
)
# Generate initial polygon coordinates
ellipse_polygons <- generate_ellipse_polygons(ellipse_params)

world <- map_data("world", wrap = c(0, 360))
world2 <- map_data("world", wrap = c(-100, 260))

worldplot <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = "gray80", color = "gray80") + 
  geom_polygon(data = world2, aes(x=long, y = lat, group = group), 
               fill = "gray80", color = "gray80") + 
  geom_polygon(
    data = ellipse_polygons,
    aes(x = Longitude, y = Latitude, group = Broad, fill = Broad, color = NA),
    alpha = 0.20,
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  # --- Pointer line: Australia Resident label to ellipse ---
  annotate(
    "segment", 
    x = 105, y = -26,        # Right edge of "Australia Resident" text
    xend = 117.5, yend = -19.3, # Outer edge of Australia ellipse
    color = "black",
    alpha = 0.40,
    linewidth = 0.4
  ) +
  # --- Pointer line: Gabon ---
  annotate(
    "segment", 
    x = 10.6, y = -3.5,        # Gabon sample
    xend = 20, yend = 3, # Gabon label
    color = "black",
    linewidth = 0.3
  ) +
  # --- Pointer line: Oman ---
  annotate(
    "segment", 
    x = 41, y = 25,        # Oman label
    xend = 55.3, yend = 17.3, # Oman sample
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(data = samp.dat, 
             aes(x = Longitude, y = Latitude, fill = Broad, shape = analysis_set), 
             size = 2, color = "black", alpha = 0.7) +
  # Stratum labels
  geom_text(data = label_data, 
            aes(x = Mean_Long, y = Mean_Lat, label = Label_Name),
            color = "black", fontface = "bold", size = 3) +
  # Placename labels
  geom_text(data = placenames, 
            aes(x = Mean_Long, y = Mean_Lat, label = Label_Name, hjust = Hjust),
            color = "black", fontface = "italic", size = 2) +
  scale_x_continuous(name = "Longitude", breaks = seq(-100, 250, by = 50), labels = lon_formatter) +
  scale_y_continuous(name = "Latitude", breaks = seq(-60, 60, by = 20), labels = lat_formatter) +
  scale_fill_manual(values = broad_stratum_colors, guide = "none") +
  scale_shape_manual(name = "Data Type", values = c("CR" = 21, "mito" = 24),
                     labels = c("CR" = "Control Region", "mito" = "Mitogenome"),
                     guide = "none") +
  coord_cartesian(xlim = c(-100, 290), ylim = c(-60, 60)) +
  theme_minimal() +
  theme(
    # 1. White background and a solid black border
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # 2. Grid lines (Thin major lines only)
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    
    # # 3. Legend placement and margins
    # legend.position = "bottom",
    # legend.box = "horizontal",
    # legend.margin = margin(t = -5, b = 0), 
    # legend.box.margin = margin(t = 0, b = 10)
  )

# --- HAWAII MAP ---
hi_inset_dat <- samp.dat |>
  filter(
    Longitude >= 189 & Longitude <= 206 &
      Latitude >= 18 & Latitude <= 26
  )

hi_labels <- tribble(
  ~Label, ~Longitude, ~Latitude,
  "NWHI", 195.5,      25.2,
  "MHI",  203.2,      22.8
)

hi_plot <- ggplot() +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), 
               fill = "gray80", color = "gray70", linewidth = 0.2) + 
  
  # REMOVED: show.legend = FALSE (allows legend generation)
  geom_point(data = hi_inset_dat, 
             aes(x = Longitude, y = Latitude, fill = Broad, shape = analysis_set), 
             size = 2, color = "black", alpha = 0.8) +
  
  geom_text(data = hi_labels, 
            aes(x = Longitude, y = Latitude, label = Label),
            fontface = "bold", size = 3, color = "black") +
  scale_x_continuous(name = "Longitude", breaks = seq(190, 205, by = 5), labels = lon_formatter) +
  scale_y_continuous(name = "Latitude", breaks = seq(18, 26, by = 4), labels = lat_formatter) +
  scale_fill_manual(values = broad_stratum_colors, guide = "none") +
  
  # Shape scale with crisp white fill override for legend keys
  scale_shape_manual(
    name = "Data Type", 
    values = c("CR" = 21, "mito" = 24),
    labels = c("CR" = "Control Region", "mito" = "Mitogenome"),
    guide = guide_legend(override.aes = list(fill = "white"))
  ) +
  
  coord_cartesian(xlim = c(189, 206), ylim = c(18, 26), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    
    # Modern ggplot2 (v3.5+) inside legend positioning syntax
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.95), # x = 2% from left, y = 95% from bottom
    legend.justification = c(0, 1),        # Anchor top-left corner of legend box
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.margin = margin(3, 5, 3, 5),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

# --- COMBINE TWO-PANEL STACKED PLOT ---
final_map <- worldplot / hi_plot + 
  plot_layout(heights = c(1, 1))

# Save PDF
pdf(file = "results-raw/mito_ocean_basin_map_ellipses.pdf", height = 6, width = 7)
final_map
dev.off()

# NAtl_samples <- filter(samp.dat, Broad == "N_Atlantic") |> 
#   mutate(weird_v_normie = ifelse(CR.hap == 50, "normie", "weird"))
# 
# world_Atl <- map_data("world")
# 
# Atl_plot <- ggplot() +
#   geom_polygon(data = world_Atl, aes(x=long, y = lat, group = group), 
#                fill = "gray60", color = "gray60") + 
#   geom_point(data = NAtl_samples, 
#              aes(x=Longitude, y = Latitude, fill = as.factor(weird_v_normie)), 
#              shape = 21, size = 2) +
#   coord_cartesian(xlim = c(-100, 0), ylim = c(0, 50)) +
#   theme_minimal() +
#   theme(panel.background = element_rect(fill = "white", color = NA),
#         panel.grid = element_line(color = "gray20"))
# 
# Atl_plot  
