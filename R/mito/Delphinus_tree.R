library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggtext)
library(tidyverse)

tree_folder <- 'aln6-all'
tree_name <- 'aln6-all_combined'
mcc_tree_file <- paste0('BEAST/xml/', tree_folder, '/', tree_name, '-MCC.tree')

mcc_tree <- treeio::read.beast(mcc_tree_file)

# Scale MCC tree branch lengths to years
mcc_tree@phylo$edge.length <- mcc_tree@phylo$edge.length * 1000000

# ==============================================================================
# Update tip labels to Scientific Name
# ==============================================================================
delphinid_map <- c(
  "EU557092_Tadu"        = "italic('Tursiops aduncus')",
  "JF289171_Fatt"        = "italic('Feresa attenuata')",
  "JN632624_Chea"        = "italic('Cephalorhynchus heavisidii')",
  "KF418381_Oorc"        = "italic('Orcinus orca')",
  "KX857264_Satt"        = "italic('Stenella attenuata')",
  "KX857383_Slon"        = "italic('Stenella longirostris')",
  "MG599457_Lhos"        = "italic('Lagenodelphis hosei')",
  "MH000365_Ddel"        = "italic('Delphinus delphis')",
  "MH992390_Lobl"        = "italic('Lagenorhynchus obliquidens')",
  "MK790767_Gmac_Shiho" = "italic('Globicephala macrorhynchus') ~ '(Shiho)'",
  "MK790780_Gmac_Naisa" = "italic('Globicephala macrorhynchus') ~ '(Naisa)'",
  "MT410901_Lalb"        = "italic('Lagenorhynchus albirostris')",
  "MT410956_Scoe"        = "italic('Stenella coeruleoalba')",
  "NC_012057_Schi"       = "italic('Sousa chinensis')",
  "NC_012059_Ttru"       = "italic('Tursiops truncatus')",
  "NC_012062_Ggri"       = "italic('Grampus griseus')",
  "NC_019590_Obre"       = "italic('Orcaella brevirostris')",
  "NC_042761_Sbre"       = "italic('Steno bredanensis')",
  "NC_045404_Steu"       = "italic('Sousa teuszii')",
  "NC_050265_Lacu"       = "italic('Lagenorhynchus acutus')",
  "NC_060611_Scly"       = "italic('Stenella clymene')",
  "NC_060612_Sfro"       = "italic('Stenella frontalis')",
  "NC_083222_Gmel"       = "italic('Globicephala melas')",
  "ON652887_Pele"        = "italic('Peponocephala electra')",
  "PP829132_Ohei"        = "italic('Orcaella heinsohni')",
  "Pcra.mito.66"        = "Pcra.mito.66",
  "Pcra.mito.55"        = "Pcra.mito.55"
)

updated_tips <- delphinid_map[mcc_tree@phylo$tip.label]
mcc_tree@phylo$tip.label <- ifelse(is.na(updated_tips), mcc_tree@phylo$tip.label, updated_tips)

# Verify updated labels
mcc_tree@phylo$tip.label

# ==============================================================================
# Generate tree
# ==============================================================================
delph_tree <- ggtree(mcc_tree) %>% revts() +
  
  # Tip labels (haplotypes)
  geom_tiplab(
    parse = TRUE,
    size = 4, 
    offset = 1000
  ) +
  
  # Label nodes with posterior support
  geom_nodelab(
    aes(
      label = case_when(
        is.na(posterior) ~ "",
        posterior >= 0.995 ~ "1",
        TRUE ~ sprintf("%.2f", posterior) # Formats < 1 to two decimal places (e.g., 0.95, 0.80)
      )
    ),    vjust = -0.5, 
    hjust = 1,
    nudge_x = -150000,
    size = 3.5, 
    color = "darkred"
  ) +
  scale_x_continuous(
    breaks = seq(-10000000, 0, 2500000), 
    minor_breaks = seq(-10000000, 0, 1250000),
    labels = function(x) format(abs(x), scientific = FALSE),
    limits = c(-10000000, 3100000) # Expands the right subpanel space past present-day (0)
  ) +
  labs(x = "Years Ago") +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10),      # Show x-axis text
        axis.ticks.x = element_line()) 

  
delph_tree

ggsave(
  filename = "results-raw/Supplemental_Delphinidae_tree.pdf", 
  plot = delph_tree, 
  width = 10, 
  height = 9, 
  onefile = FALSE
)
