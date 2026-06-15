library(tidyverse)
library(ggplot2)
library(ggrepel)
library(maps)
library(mapdata)
library(sp)
library(scales)

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
         Broad = ifelse(Broad %in% c('MHI', 'NWHI'), 'HI_res', Broad),
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
                          Broad %in% c("NPac_Eastern", "NPac_Central", "S_Pacific", "HI_res"), 
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
  "HI_res"
)

# Update samp.dat to make 'Broad' a factor with those levels
samp.dat <- samp.dat |>
  # this filters out two Korean market samples and a HI sample that could be
  # MHI resident or NWHI resident
  filter(Broad != 'Unknown') |> 
  # Move HI_res samples to the end of the data frame so they plot last
  arrange(Broad == "HI_res") |> 
  arrange(analysis_set == "mito") |> 
  mutate(Broad = factor(Broad, levels = basin_order))

# Create a named vector for the "pretty" labels
basin_labels <- c(
  "N_Atlantic"    = "North Atlantic",
  "S_Atlantic"    = "South Atlantic",
  "S_Pacific"     = "South Pacific",
  "NPac_Western"  = "Western North Pacific",
  "NPac_Eastern"  = "Eastern North Pacific",
  "NPac_Central"  = "Central North Pacific",
  "Indian_Ocean"  = "Indian Ocean",
  "Aus_res"       = "Australia Resident",
  "HI_res"        = "Hawaii Resident"
)

# Calculate the center point for each stratum to place the label
# label_data <- samp.dat |>
#   group_by(Broad) |>
#   summarize(
#     Mean_Long = mean(Longitude, na.rm = TRUE),
#     Mean_Lat = mean(Latitude, na.rm = TRUE)
#   ) |>
#   mutate(Label_Name = basin_labels[as.character(Broad)])
label_data <- tribble(
  ~Broad,          ~Mean_Long,  ~Mean_Lat,  ~Label_Name,
  "N_Atlantic",    -50,    25.14908,  "North\nAtlantic",
  "S_Atlantic",    -10,   -20,  "South\nAtlantic",
  "Indian_Ocean",   85.37323,     -10,  "Indian Ocean",
  "S_Pacific",     230,   -25,  "South\nPacific",
  "NPac_Western",  144.49720,    15.75972,  "Western\nN. Pacific",
  "NPac_Central",  197.05908,    40,  "Central\nN. Pacific",
  "NPac_Eastern",  264.44204,    14.89037,  "Eastern\nN. Pacific",
  "Aus_resident",  124.41908,   -15.90395,  "Australia Resident",
  "HI_res",        202.30489,    20.74375,  "MHI & NWHI\nResidents"
)
world <- map_data("world", wrap = c(0, 360))
world2 <- map_data("world", wrap = c(-100, 260))

worldplot <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = "gray60", color = "gray60") + 
  geom_polygon(data = world2, aes(x=long, y = lat, group = group), 
               fill = "gray60", color = "gray60") + 
  geom_point(data = samp.dat, 
             aes(x = Longitude, y = Latitude, fill = Broad, shape = analysis_set), 
             size = 2, color = "black", alpha = 0.7) +
  geom_text(data = label_data, 
            aes(x = Mean_Long, y = Mean_Lat, label = Label_Name),
            color = "black", fontface = "bold", size = 3) +
  scale_x_continuous(name = "Longitude", breaks = seq(-100, 250, by = 50), labels = lon_formatter) +
  scale_y_continuous(name = "Latitude", breaks = seq(-60, 60, by = 20), labels = lat_formatter) +
  scale_fill_manual(values = broad_stratum_colors, guide = "none") +
  scale_shape_manual(name = "Data Type", values = c("CR" = 21, "mito" = 24),
                     labels = c("CR" = "Control Region", "mito" = "Mitogenome")) +
  coord_cartesian(xlim = c(-100, 290)) +
  theme_minimal() +
  theme(
    # 1. White background and a solid black border
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # 2. Grid lines (Thin major lines only)
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    
    # 3. Legend placement and margins
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = -5, b = 0), 
    legend.box.margin = margin(t = 0, b = 10)
  )
pdf(file = "results-raw/mito_ocean_basin_map.pdf", height = 6, width = 7)
worldplot
dev.off()

NAtl_samples <- filter(samp.dat, Broad == "N_Atlantic") |> 
  mutate(weird_v_normie = ifelse(CR.hap == 50, "normie", "weird"))

world_Atl <- map_data("world")

Atl_plot <- ggplot() +
  geom_polygon(data = world_Atl, aes(x=long, y = lat, group = group), 
               fill = "gray60", color = "gray60") + 
  geom_point(data = NAtl_samples, 
             aes(x=Longitude, y = Latitude, fill = as.factor(weird_v_normie)), 
             shape = 21, size = 2) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(0, 50)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "gray20"))

Atl_plot  
