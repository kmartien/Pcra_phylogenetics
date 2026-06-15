library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggpattern)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('../Pcra.database.data/data/ASPhyloCR.rda')
load('data/broad_stratum_colors.rda')

CR.haps <- filter(CR.haps, !is.na(CR.hap))

#tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc_combined/aln.7.CDS-rRNA_constant_orc_combined-MCC.trees')
#tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc/aln.3.CDS-rRNA_constant_orc-MCC.tree')
tree <- read.beast(file = 'BEAST/xml/Pcra_CR_unique_constant_orc_logMRCA/Pcra_CR_unique_constant_orc_logMRCA_combined_MCC.tree')

CR.strata <- 
  filter(CR.haps, Animal.ID %in% ASPhyloCR$Animal.ID) |> 
  left_join(mito.haps) |> 
  left_join(Pcra.strata) |> 
  mutate(Broad = ifelse(Broad %in% c('Atl-weird', 'Atl-normal', 'Atl-unknown'), 'N_Atlantic', Broad),
         Broad = ifelse(Broad %in% c('MHI', 'NWHI'), 'HI_res', Broad),
         Broad = ifelse(Broad == 'Atl-South', 'S_Atlantic', Broad),
         Broad = ifelse(Broad == 'Spac', 'S_Pacific', Broad),
         Broad = ifelse(Broad == 'Wpac', 'NPac_Western', Broad),
         Broad = ifelse(Broad == 'Epac', 'NPac_Eastern', Broad),
         Broad = ifelse(Broad == 'CNP', 'NPac_Central', Broad),
         Broad = ifelse(Broad == 'Indian Ocean', 'Indian_Ocean', Broad),
         # Change stratum for sample that could be either MHI or NWHI
         Broad = ifelse(Animal.ID == 'z0203251', 'HI_res', Broad))

tree@phylo$edge.length <- tree@phylo$edge.length*1000000

##############################################################################
# Color-code tips by broad strata

Broad_strata <- CR.strata |> 
  select(CR.hap, Broad) |> 
  distinct() |> 
  rename(label = CR.hap) |> 
  # Create an x-position for each dot within a haplotype
  group_by(label) |> 
  mutate(dot_x = seq_along(Broad)) |> 
  ungroup() |> 
  # Fix factor levels to include HI_res and Aus_resident
  mutate(Broad = factor(Broad, levels = c(
    "N_Atlantic", "S_Atlantic", "Indian_Ocean", "S_Pacific",
    "NPac_Western", "NPac_Central", "NPac_Eastern", 
    "Aus_resident", "HI_res", "Unknown"
  )))

p <- ggtree(tree, mrsd = '2025-01-01') +
  geom_tree() + theme_tree2() + 
#  geom_text(aes(label=node)) +
  geom_tiplab() +
  scale_linetype_manual(values = c(2,1)) +
  geom_point2(aes(subset = posterior > 0.8), size = 2) +
  scale_x_continuous(breaks = seq(-225000, 0, 50000), minor_breaks = seq(-225000, 0, 10000)) +
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
p <- p +
  geom_facet(panel = "Ocean Basin", 
             data = Broad_strata, 
             geom = geom_point,
             # Using 0.8 and 1.2 forces the dots to be much closer together
             mapping = aes(x = ifelse(dot_x == 1, 0.8, 1.2), fill = Broad), 
             shape = 21,   
             size = 3,   # Increased slightly so they touch
             color = "black",
             stroke = 0.5) + 
  scale_fill_manual(values = broad_stratum_colors, name = "Basin") +
  
  # FIX: Use xlim_expand for faceted panels instead of xlim_tree
  xlim_expand(c(0, 2), panel = "Ocean Basin") + 
  
  theme(
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines"), 
    plot.margin = margin(10, 20, 10, 10),
    legend.position = "right",
    text = element_text(size = 14)
  )

# Keep a very small width for the dot panel to force compression
p <- facet_widths(p, c(9.5, 0.5))

pdf(file = 'results-raw/CR_viridis_tree_logMRCA.pdf', width = 8, height = 7)
p
dev.off()

