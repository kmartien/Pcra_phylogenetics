library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggbreak)
library(tidyverse)

tree_folder <- 'aln3-all'
tree_name <- 'aln3-all_skyline_combined'
tree_posterior_file <- paste0('BEAST/xml/', tree_folder, '/', tree_name, '-Tree.trees')
mcc_tree_file <- paste0('BEAST/xml/', tree_folder, '/', tree_name, '-MCC.tree')
mcc_tree_layout <- 'rectangular'
densitree_layout <- 'slanted'
min_posterior_support <- 0.94

# ==============================================================================
# 1. LOAD DATA AND PREPROCESS STRATA
# ==============================================================================
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

mito.strata <- left_join(mito.haps, CR.haps) |> 
  left_join(Pcra.strata) |> 
  mutate(Broad = ifelse(Broad %in% c('Atl-weird', 'Atl-normal', 'Atl-unknown'), 'N_Atlantic', Broad),
#         Broad = ifelse(Broad %in% c('MHI', 'NWHI'), 'HI_res', Broad),
         Broad = ifelse(Broad == 'Atl-South', 'S_Atlantic', Broad),
         Broad = ifelse(Broad == 'Spac', 'S_Pacific', Broad),
         Broad = ifelse(Broad == 'Wpac', 'NPac_Western', Broad),
         Broad = ifelse(Broad == 'Epac', 'NPac_Eastern', Broad),
         Broad = ifelse(Broad == 'CNP', 'NPac_Central', Broad),
         Broad = ifelse(Broad == 'Indian Ocean', 'Indian_Ocean', Broad))

Broad_strata <- mito.strata |> 
  select(mito.hap, Broad) |> 
  distinct() |> 
  rename(label = mito.hap) |> 
  group_by(label) |> 
  mutate(dot_x = seq_along(Broad)) |> 
  ungroup() |> 
  mutate(Broad = factor(Broad, levels = c(
    "N_Atlantic", "S_Atlantic", "Indian_Ocean", "S_Pacific",
    "NPac_Western", "NPac_Central", "NPac_Eastern", 
    "Aus_resident", "MHI", "NWHI", "Unknown"
  ))) |> 
  mutate(dot_x = as.numeric(dot_x))

strata_labels <- tribble(
  ~Broad, ~Name,
  "N_Atlantic", "North Atlantic",
  "S_Atlantic", "South Atlantic",
  "Indian_Ocean", "Indian Ocean",
  "S_Pacific", "South Pacific",
  "NPac_Western", "Western N. Pacific",
  "NPac_Central", "Central N. Pacific",
  "NPac_Eastern", "Eastern N. Pacific",
  "Aus_resident", "Australia Resident",
  "MHI", "MHI Resident",
  "NWHI", "NWHI Resident",
  "Unknown", "Unknown"
)

offset_val <- 5000 # Horizontal spacing for tip dots

# ==============================================================================
# 2. READ & PROCESS POSTERIOR TREES (DENSITREE CLOUD)
# ==============================================================================
raw_posterior <- ape::read.nexus(tree_posterior_file)

burn_in_fraction <- 0.10
total_trees <- length(raw_posterior)
burn_in_cutoff <- round(total_trees * burn_in_fraction)
post_burn_in_trees <- raw_posterior[(burn_in_cutoff + 1):total_trees]

subsample_size <- 100
sampled_indices <- round(seq(1, length(post_burn_in_trees), length.out = subsample_size))
final_trees <- post_burn_in_trees[sampled_indices]

# Scale posterior tree branch lengths to years
final_trees <- lapply(final_trees, function(t) {
  t$edge.length <- t$edge.length * 1000000
  return(t)
})
class(final_trees) <- "multiPhylo"

# ==============================================================================
# 3. READ & PROCESS THE MCC SUMMARY TREE
# ==============================================================================
mcc_tree <- treeio::read.beast(mcc_tree_file)

# # Cladewise edge-cascade to collapse unsupported nodes (< min_posterior_support) into polytomies
# mcc_tree@phylo <- ape::reorder.phylo(mcc_tree@phylo, "cladewise")
# low_support_nodes <- mcc_tree@data %>%
#   filter(!is.na(posterior) & posterior < min_posterior_support) %>%
#   pull(node)
# 
# for (i in 1:nrow(mcc_tree@phylo$edge)) {
#   parent_node <- mcc_tree@phylo$edge[i, 1]
#   child_node  <- mcc_tree@phylo$edge[i, 2]
#   if (child_node %in% low_support_nodes) {
#     len_to_push <- mcc_tree@phylo$edge.length[i]
#     if (len_to_push > 0) {
#       child_edge_indices <- which(mcc_tree@phylo$edge[, 1] == child_node)
#       mcc_tree@phylo$edge.length[child_edge_indices] <- mcc_tree@phylo$edge.length[child_edge_indices] + len_to_push
#       mcc_tree@phylo$edge.length[i] <- 0
#     }
#   }
# }

# Scale MCC tree branch lengths to years
mcc_tree@phylo$edge.length <- mcc_tree@phylo$edge.length * 1000000

# Fortify MCC tree layout into data frame coordinates
mcc_df <- ggtree::fortify(mcc_tree, layout = mcc_tree_layout)
mcc_tip_order <- mcc_df %>% filter(isTip) %>% arrange(y) %>% pull(label)

# Shift MCC coordinates leftward so modern tips sit exactly at 0
mcc_max_x <- max(mcc_df$x, na.rm = TRUE)
mcc_df$x <- mcc_df$x - mcc_max_x

# Map shifted coordinates to tip data frame for basin color dots
tip_coords <- mcc_df %>% 
  filter(isTip) %>% 
  select(label, x, y)

plot_data <- Broad_strata %>%
  inner_join(tip_coords, by = "label")

# ==============================================================================
# 4. DEFINE CLIMATE CHRONOLOGY TIMELINE (MIS)
# Based on The LR04 Benthic Stack (Lisiecki & Raymo, 2005)
# The Optimized MIS Chronology (Railsback et al., 2015) is a new authoritative
# source, but it just subdivides MISs from Lisiecki & Raymo into smaller sub-
# stages, which I'm not using
# ==============================================================================
mis_timeline <- data.frame(
  xmin = c(-14000, -29000, -57000, -71000, -130000, -191000, -243000, -300000, -337000, -374000, -424000, -478000, -524000),
  xmax = c(0,      -14000, -29000, -57000, -71000,  -130000, -191000, -243000, -300000, -337000, -374000, -424000, -478000),
  stage = c("", "MIS\n2", "MIS\n3", "MIS\n4", "MIS\n5", "MIS\n6", "MIS\n7", "MIS\n8", "MIS\n9", "MIS\n10", "MIS\n11", "MIS\n12", "MIS\n13"),
  type = c("interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial", 
           "interglacial", "glacial", "interglacial")
) %>% 
  mutate(x_center = (xmin + xmax) / 2)

# ==============================================================================
# 5. INITIALIZE PLOT AND DEFINE LAYERS
# ==============================================================================
# Initialize base cloud synchronized to the MCC tree's tip positions
p <- ggdensitree(
  final_trees, 
  alpha = 0.1, 
  colour = "steelblue", 
  layout = densitree_layout, #  rectangular or slanted 
  tip.order = mcc_tip_order
  )

# Layer 1: Shaded Glacial Bands (Absolute Bottom)
rect_layer <- geom_rect(
  data = dplyr::filter(mis_timeline, type == "glacial"),
  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  fill = "grey93",
  alpha = 1.0,
  inherit.aes = FALSE
)

# Manual Grid Line Layers (Sandwiched right above the gray rectangles)
minor_breaks_vec <- seq(-450000, 0, 10000)
major_breaks_vec <- seq(-450000, 0, 50000)

manual_minor_grid <- geom_vline(
  xintercept = minor_breaks_vec, 
  color = "grey80", 
  linewidth = 0.2, 
  inherit.aes = FALSE
)

manual_major_grid <- geom_vline(
  xintercept = major_breaks_vec, 
  color = "grey75", # Slightly darker than minor grid to stand out cleanly
  linewidth = 0.3, 
  inherit.aes = FALSE
)

# Layer 3: The Consensus MCC Tree (Solid Black)
mcc_layer <- geom_tree(
  data = mcc_df,
  colour = "black",
  size = 0.8,
  layout = "rectangular",
  inherit.aes = FALSE
)

# Layer 4: Strongly Supported Node Dots (> min_posterior_support) on the MCC tree
node_points_layer <- geom_point(
  data = dplyr::filter(mcc_df, !isTip & posterior > min_posterior_support),
  aes(x = x, y = y),
  shape = 21,
  fill = "black",
  color = "black",
  size = 2,
  inherit.aes = FALSE
)

# Layer 5: Color-Coded Tip Dots by Ocean Basin
tip_points_layer <- geom_point(
  data = plot_data, 
  aes(x = x + (dot_x * offset_val), y = y, fill = Broad),
  shape = 21, 
  size = 2.5, 
  color = "black", 
  stroke = 0.5,
  inherit.aes = FALSE 
)

# Layer 6: Climate Timeline Labels (Top-most layer)
text_layer <- geom_text(
  data = mis_timeline,
  aes(x = x_center, y = Inf, label = stage),
  vjust = 2,
  size = 3,
  fontface = "bold",
  color = "grey40",
  inherit.aes = FALSE
)

# ==============================================================================
# 6. ASSEMBLE LAMINATED CANVAS & APPLY FINAL COSMETICS
# ==============================================================================
# Order components sequentially from back to front:
# [Rectangles] -> [Minor Grid] -> [Major Grid] -> [Densitree Cloud] -> [MCC Tree] -> [Nodes] -> [Tips] -> [Labels]
p$layers <- c(rect_layer, manual_minor_grid, manual_major_grid, p$layers, mcc_layer, node_points_layer, tip_points_layer, text_layer)

# Apply scales, viewport crop, and clean theme
p <- p +
  scale_fill_manual(values = broad_stratum_colors, labels = deframe(strata_labels), name = "Basin") +
  theme_tree2() +
  scale_x_continuous(
    breaks = major_breaks_vec, 
    minor_breaks = minor_breaks_vec,
    labels = function(x) format(round(x), scientific = FALSE),
    expand = expansion(mult = c(0.0, 0.05)) # Left expansion set to 0 to snap to -300,000 limit
  ) +
  # Non-destructive viewport crop (Truncates timeline at -300,000 without deleting tree roots)
  coord_cartesian(xlim = c(-465000, 1000)) +
  theme(
    # Completely turn off native theme grid lines since we drew our own manual sandwich grid
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(size = 10),
    axis.ticks.x       = element_line(),
    legend.position    = c(0.15, 0.60), # x = 12% from left margin, y = 70% up from bottom
    legend.background  = element_rect(colour = "black"), # No ugly white box outline
    legend.key         = element_rect(fill = "transparent", colour = NA), # No grey boxes under legend icons
    legend.key.height  = unit(10, "pt"),  # Shrinks the vertical bounding box of the colored dot
    legend.spacing.y   = unit(3, "pt"),   # Sets a tight gap between adjacent rows
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 9)
  ) +
  labs(
    x = "Years Ago"
  )

# ==============================================================================
# 7. EXPORT TO FILE
# ==============================================================================
pdf(file = paste0('results-raw/', tree_name, '-mcc+densitree_', mcc_tree_layout, '-', densitree_layout, '_temp.pdf'), width = 8, height = 6)
print(p)
dev.off()

# ==============================================================================
# 8. TREE WITH LABELED NODES FOR SUPPLEMENTAL MATERIAL
# ==============================================================================
# Define exactly which ticks you want visible on both sides of the break
my_breaks <- c(-380000, -150000, -100000, -50000, 0)

suppl_p <- ggtree(mcc_tree) %>% revts() +
  # Break axis between -350kya and -150kya
  scale_x_break(c(-350000, -150000), symbol = "slash") +
  
  # Tip labels (haplotypes)
  geom_tiplab(size = 3, offset = 1000) + 
  
  # Internal node numbers
  geom_label2(
    aes(label = node, subset = !isTip),
    fill = "white",          # White box masks the branch line
    color = "darkred",       # Text color
    label.size = 0,          # Removes the border box outline
    label.padding = unit(0.1, "lines"), # Keeps the white box tight
    size = 2.5,
    vjust = 0.5,             # Centered vertically on the node
    hjust = 0.5              # Centered horizontally on the node
  ) +
  
  # Custom tick positions + right coordinate space reserved for tip text
  scale_x_continuous(
    breaks = my_breaks,
    labels = function(x) format(abs(x), scientific = FALSE),
    limits = c(-390000, 15000) # Expands the right subpanel space past present-day (0)
  ) +
  
  theme_tree2() +
  labs(x = "Years Ago") +
  
  # SUPPRESS TOP X-AXIS: Remove the duplicate top axis generated by ggbreak
  theme(
    axis.text.x.top  = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank()
  )

# ==============================================================================
# EXPORT PDF (onefile = FALSE kills the extra blank page 1)
# ==============================================================================
ggsave(
  filename = "results-raw/Supplemental_MCC_Tree_Break.pdf", 
  plot = suppl_p, 
  width = 9, 
  height = 10, 
  onefile = FALSE
)
# For ggbreak, cite:
# S Xu#, M Chen#, T Feng, L Zhan, L Zhou, G Yu*. Use ggbreak to effectively utilize plotting space to deal with large datasets and outliers. Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846
# Supplemental file: https://github.com/YuLab-SMU/supplemental-ggbreak.

# ==============================================================================
# 9. GENERATE SUPPLEMENTAL NODE TABLE
# ==============================================================================

# Convert tree object to a data frame
mcc_tbl <- as_tibble(mcc_tree)

# Get the total number of tip nodes (internal nodes are > num_tips)
num_tips <- Ntip(mcc_tree@phylo)

# Dynamically identify the HPD column name generated by TreeAnnotator
hpd_col_name <- grep("height_.*HPD", names(mcc_tbl), value = TRUE)[1]

suppl_node_table <- mcc_tbl %>%
  filter(node > num_tips) %>%  # Isolate internal nodes only
  rowwise() %>%
  mutate(
    # Scale divergence times to years (matching your 1e6 tree scaling factor)
    divergence_time_years = height * 1e6,
    
    # Safely extract lower and upper 95% HPD bounds from the list-column
    hpd_vec = list(if (!is.null(hpd_col_name)) get(hpd_col_name) else c(NA, NA)),
    lower_95_HPD = hpd_vec[1] * 1e6,
    upper_95_HPD = hpd_vec[2] * 1e6
  ) %>%
  ungroup() %>%
  select(
    node,
    posterior,
    divergence_time_years,
    lower_95_HPD,
    upper_95_HPD
  ) %>%
  mutate(
    posterior = round(posterior, 3),
    divergence_time_years = round(divergence_time_years, 0),
    lower_95_HPD = round(lower_95_HPD, 0),
    upper_95_HPD = round(upper_95_HPD, 0)
  ) %>%
  arrange(node)

# Save summary table as a CSV
write_csv(suppl_node_table, "results-raw/Supplemental_Table_Node_Support.csv")

