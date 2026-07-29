library(ape)
library(ggtree)
library(ggplot2)
library(tidyverse)

# ==============================================================================
# 1. DEFINE CLEAN METADATA FOR YOUR 7 NUCLEAR GENOMES
# ==============================================================================
sample_metadata <- data.frame(
  label = c('Oorc_Morgan', 
            'Pcra_z0027510', 'Pcra_z0018462', 'Pcra_z0045928', 
            'Gmac_z0017981', 'Gmac_z0112653', 'Gmac_z0112731'),
  clean_name = c('Orca Outgroup', 
                 'FKW N Atl', 'FKW ETP', 'FKW MHI', 
                 'SFPW ETP', 'SFPW MHI', 'SFPW N Atl')
)

# ==============================================================================
# 2. READ THE REARRANGED NUCLEAR BUSCO TREE
# ==============================================================================
nuc_tree <- ape::read.tree('hPSMC_Dsuite/Dsuite/Busco_phylogeny/Pcra_Gmac_Oorc_129_largest_Busco_rearranged.newick')

# ==============================================================================
# 3. BUILD THE BASE TREE PANEL (AXIS REMOVED)
# ==============================================================================
p_tree <- ggtree(nuc_tree) %<+% sample_metadata +
  geom_tree(linewidth = 0.8) +
  # Render clean geographical names right next to the tree tips
  geom_tiplab(aes(label = clean_name), offset = 0.0001, size = 3.5, fontface = "bold") +
  theme_tree()

# Determine the tree's maximum length so our offsets scale perfectly
max_tree_x <- max(p_tree$data$x, na.rm = TRUE)

# ==============================================================================
# 4. PARSE FBRANCH DATA INTO GHEATMAP FORMAT (ALL 6 RECIPIENTS)
# ==============================================================================
fb_raw <- read.delim('results-raw/Dsuite/Pcra_dsuite_results_combined_Fbranch_noZscore.txt', sep='\t')

fb_matrix_df <- fb_raw %>%
  filter(branch_descendants %in% nuc_tree$tip.label) %>% # Verified %in% syntax
  column_to_rownames("branch_descendants") %>%
  # Select all 6 target recipient lineages (3 FKW and 3 SFPW)
  select(Pcra_z0027510, Pcra_z0018462, Pcra_z0045928, 
         Gmac_z0017981, Gmac_z0112653, Gmac_z0112731) %>%
  # Rename to clean headers that differentiate species and basin
  rename(`FKW N Atl` = Pcra_z0027510, 
         `FKW ETP`   = Pcra_z0018462, 
         `FKW MHI`   = Pcra_z0045928,
         `SFPW ETP`  = Gmac_z0017981,
         `SFPW MHI`  = Gmac_z0112653,
         `SFPW N Atl`  = Gmac_z0112731)

# ==============================================================================
# 5. LAMINATE THE WIDE HEATMAP & STRIP WHITESPACE
# ==============================================================================
p_final <- gheatmap(
  p = p_tree, 
  data = fb_matrix_df, 
  offset = max_tree_x * 0.22,  
  width = 1.6,                 
  color = "gray70",             
  colnames_position = "top",   
  colnames_angle = 45,         
  colnames_offset_y = -0.25,   # Negative value counteracts default padding to pull labels tight
  hjust = 0                    
) +
  # scale_fill_viridis_c(
  #   option = "rocket", 
  #   direction = -1, 
  #   na.value = "grey95", 
  #   name = expression(italic(f)[branch])
  # ) +
  # Custom gradient injection: Forces 0 to be pure white, then seamlessly blends into rocket
  scale_fill_gradientn(
    colors = c("white", scales::viridis_pal(option = "rocket", direction = -1)(100)),
    limits = c(0, NA),
    na.value = "grey95", 
    name = expression(italic(f)[branch])
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  # Expands the native Y-axis space up by 30%, holding the 45-degree text inside the layout window
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  # Allows the rightmost column name to project past the grid boundary without being clipped
  coord_cartesian(clip = "off") + 
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.1),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    # USER OPTIMIZED: Your perfect tight legend spacing
    legend.box.spacing = unit(5, "pt"), 
    # Standard clean padding
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

# ==============================================================================
# 6. EXPORT GRAPHIC
# ==============================================================================
pdf(file = 'results-raw/nuclear_fbranch_1to1_perfect_sort.pdf', width = 8.5, height = 5)
print(p_final)
dev.off()