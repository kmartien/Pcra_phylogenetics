library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)

tree_name <- 'aln3-all'
tree_posterior_file <- paste0('BEAST/xml/', tree_name, '/', tree_name, '_combined-Tree.trees')
mcc_tree_file <- paste0('BEAST/xml/', tree_name, '/', tree_name, '_combined-MCC.tree')
mcc_tree_layout <- 'rectangular'
densitree_layout <- 'slanted'
min_posterior_support <- 0.8

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
         Broad = ifelse(Broad %in% c('MHI', 'NWHI'), 'HI_res', Broad),
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
    "Aus_resident", "HI_res", "Unknown"
  ))) |> 
  mutate(dot_x = as.numeric(dot_x))

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
# ==============================================================================
mis_timeline <- data.frame(
  xmin = c(-14000, -29000, -57000, -71000, -130000, -191000, -243000, -300000),
  xmax = c(0,      -14000, -29000, -57000, -71000,  -130000, -191000, -243000),
  stage = c("", "MIS\n2", "MIS\n3", "MIS\n4", "MIS\n5", "MIS\n6", "MIS\n7", "MIS\n8"),
  type = c("interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial")
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

# NEW: Manual Grid Line Layers (Sandwiched right above the gray rectangles)
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
  size = 3, 
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
  scale_fill_manual(values = broad_stratum_colors, name = "Basin") +
  theme_tree2() +
  scale_x_continuous(
    breaks = major_breaks_vec, 
    minor_breaks = minor_breaks_vec,
    labels = function(x) format(round(x), scientific = FALSE),
    expand = expansion(mult = c(0.0, 0.05)) # Left expansion set to 0 to snap to -300,000 limit
  ) +
  # NEW: Non-destructive viewport crop (Truncates timeline at -300,000 without deleting tree roots)
  coord_cartesian(xlim = c(-465000, 1000)) +
  theme(
    # Completely turn off native theme grid lines since we drew our own manual sandwich grid
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(size = 10),
    axis.ticks.x       = element_line(),
    legend.position    = c(0.30, 0.70), # x = 12% from left margin, y = 70% up from bottom
    legend.background  = element_rect(colour = "black"), # No ugly white box outline
    legend.key         = element_rect(fill = "transparent", colour = NA), # No grey boxes under legend icons
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9)
  ) +
  labs(
    x = "Years Ago"
  )

# ==============================================================================
# 7. EXPORT TO FILE
# ==============================================================================
pdf(file = paste0('results-raw/', tree_name, '-mcc+densitree_', mcc_tree_layout, '-', densitree_layout, '.pdf'), width = 8, height = 7)
print(p)
dev.off()