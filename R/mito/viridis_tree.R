library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggpattern)
library(tidyverse)
library(viridis)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

#tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc_combined/aln.7.CDS-rRNA_constant_orc_combined-MCC.trees')
#tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc/aln.3.CDS-rRNA_constant_orc-MCC.tree')
tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc_logMRCA/aln.3.CDS-rRNA_constant_orc_logMRCA-MCC.tree')

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

################################################################################
# # Collapse nodes with support < 0.8
# # 1. Reorder the tree to "cladewise" (pre-order layout). 
# # This guarantees parent rows always appear BEFORE their child rows in the edge matrix.
# tree@phylo <- ape::reorder.phylo(tree@phylo, "cladewise")
# 
# # 2. Identify internal nodes with posterior support less than 0.8
# low_support_nodes <- tree@data %>%
#   filter(!is.na(posterior) & posterior < 0.8) %>%
#   pull(node)
# 
# # 3. Loop through every single edge in the tree sequentially
# for (i in 1:nrow(tree@phylo$edge)) {
#   parent_node <- tree@phylo$edge[i, 1]
#   child_node  <- tree@phylo$edge[i, 2]
#   
#   # If the edge leads TO a low-support node, we collapse it
#   if (child_node %in% low_support_nodes) {
#     len_to_push <- tree@phylo$edge.length[i]
#     
#     if (len_to_push > 0) {
#       # Find all immediate downstream branches stemming from this child node
#       child_edge_indices <- which(tree@phylo$edge[, 1] == child_node)
#       
#       # Cascade the length down to the next generation of branches
#       tree@phylo$edge.length[child_edge_indices] <- tree@phylo$edge.length[child_edge_indices] + len_to_push
#       
#       # Set the collapsed parent branch length to exactly 0
#       tree@phylo$edge.length[i] <- 0
#     }
#   }
# }
# 
# tree@phylo$edge.length[edge_rows] <- 0

################################################################################

# Make divergence time and uncertainty be in calendar years instead of Myears
tree@phylo$edge.length <- tree@phylo$edge.length*1000000
tree@data <- tree@data %>%
  mutate(height_0.95_HPD = map2(height_0.95_HPD, posterior, function(hpd_vector, post_value) {
    # 1. Drop if no HPD data exists
    # 2. Drop if posterior is NA (which happens at the tips)
    # 3. Drop if posterior support is 0.8 or lower
    if (is.null(hpd_vector) || is.na(post_value) || post_value <= 0.8) {
      return(NULL)
    }
    # Otherwise, scale the HPD vector by 1,000,000
    return(hpd_vector * 1000000)
  }))

##############################################################################

# Color-code tips by broad strata

Broad_strata <- mito.strata |> 
  select(mito.hap, Broad) |> 
  distinct() |> 
  rename(label = mito.hap) |> 
  # Create an x-position for each dot within a haplotype
  group_by(label) |> 
  mutate(dot_x = seq_along(Broad)) |> 
  ungroup() |> 
  # Fix factor levels to include HI_res and Aus_resident
  mutate(Broad = factor(Broad, levels = c(
    "N_Atlantic", "S_Atlantic", "Indian_Ocean", "S_Pacific",
    "NPac_Western", "NPac_Central", "NPac_Eastern", 
    "Aus_resident", "HI_res", "Unknown"
  ))) |> 
  mutate(dot_x = as.numeric(dot_x))

clade_colors <- turbo(n = 9)[c(2,4,5,6,7,9)]

offset_val <- 5000 # how far second dot should be moved for haps with multiple strata

# Defining glacial intervals using negative X-axis scale
# values from Lisiecki and Raymo (2005) and Railsback et al. (2015)
# 1. Define the complete timeline for ALL Marine Isotope Stages (MIS)
mis_timeline <- data.frame(
  xmin = c(-14000, -29000, -57000, -71000, -130000, -191000, -243000, -300000),
  xmax = c(0,      -14000, -29000, -57000, -71000,  -130000, -191000, -243000),
  stage = c("", "MIS\n2", "MIS\n3", "MIS\n4", "MIS\n5", "MIS\n6", "MIS\n7", "MIS\n8"),
  type = c("interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial", "interglacial", "glacial")
) %>% 
  # Mathematically find the center of each band for text placement
  mutate(x_center = (xmin + xmax) / 2)

p <- ggtree(tree, mrsd = '2025-01-01') %<+% Broad_strata +
  # BACKGROUND LAYER: Shaded glacial bands
  geom_rect(
    data = dplyr::filter(mis_timeline, type == "glacial"),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "grey90",      # Soft light gray
    alpha = 0.5,          # Semi-transparent
    inherit.aes = FALSE
  ) +

  # LAYER 2: Chronology labels pinned to the very top inside each band
  geom_text(
    data = mis_timeline,
    aes(x = x_center, y = Inf, label = stage),
    vjust = 1,            # Pushes the text down slightly so it sits inside the top margin
    size = 3,             # Font size
    fontface = "bold",    # Clean, scannable bold text
    color = "grey40",     # Muted dark gray so it doesn't distract from tip labels
    inherit.aes = FALSE
  ) +
  
  geom_tree() + 
  theme_tree2() 
  ## this code adds colored bars at the ends of clades - not using it for now
  # geom_cladelab(node = 99, label = "", barcolor = clade_colors[1]) +
  # geom_cladelab(node = 123, label = "", barcolor = clade_colors[2]) +
  #  geom_text(aes(label=node)) +

# Extract the coordinate mapping from the plot object
# This table contains 'label', 'x', and 'y' for every tip
tip_coords <- p$data %>% 
  filter(isTip) %>% 
  select(label, x, y)

# Join these coordinates to your Broad_strata data
# This ensures every instance of a haplotype gets the correct Y-axis position
plot_data <- Broad_strata %>%
  inner_join(tip_coords, by = "label") %>%
  mutate(dot_x = as.numeric(dot_x))

# Add the points using a standard geom_point layer
p <- p + 
  geom_point(
    data = plot_data, 
    aes(
      # We use the tip's 'x' coordinate plus the offset
      x = x + (dot_x * offset_val), 
      y = y, 
      fill = Broad
    ),
    shape = 21, 
    size = 3, 
    color = "black", 
    stroke = 0.5,
    # This prevents the layer from trying to find 'node' or other tree vars
    inherit.aes = FALSE 
  ) +
  scale_fill_manual(values = broad_stratum_colors, name = "Basin") +
  # Expand the x-axis so the dots aren't cut off
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))

p <- p +
  scale_linetype_manual(values = c(2,1)) +
  geom_point2(aes(subset = posterior > 0.8), size = 2) +
  geom_range(range = "height_0.95_HPD", color = "darkred", size = 1.5, alpha = 0.5) +
  scale_x_continuous(breaks = seq(-300000, 0, 50000), minor_breaks = seq(-300000, 0, 10000)) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10),      # Show x-axis text
        axis.ticks.x = element_line()) 

##################
# Try this for collapsing unsupported nodes into polytomies
#as.polytomy(tree, feature='node.label', fun=function(x) as.numeric(x) < 70)
##################
 
#p <- flip(p, 73, 85)

pdf(file = 'results-raw/viridis_tree_logMRCA.pdf', width = 8, height = 7)
p
dev.off()

