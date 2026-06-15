library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggpattern)
library(tidyverse)
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

tree@phylo$edge.length <- tree@phylo$edge.length*1000000

##############################################################################
# Color-code tips by broad strata

Broad_strata <- select(mito.strata, c('mito.hap', 'Broad')) |> 
  distinct() |> 
  rename(label = mito.hap) |> 
  group_by(label) |> 
  mutate(proportion = 1/n()) |> 
  ungroup() |> 
#  mutate(pattern_type = ifelse(Broad %in% c("MHI", "NWHI", "Aus_resident"), "crosshatch", "none")) |> 
  mutate(Broad = factor(Broad, levels = c("N_Atlantic", "S_Atlantic", 
                                          "Indian_Ocean", "S_Pacific",
                                          "NPac_Western", "NPac_Central",
                                          "NPac_Eastern", "HI_res",
                                          "Aus_resident", "Unknown")))

p <- ggtree(tree, mrsd = '2025-01-01') +
  geom_tree() + theme_tree2() + 
#  geom_text(aes(label=node)) +
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
as.polytomy(tree, feature='node.label', fun=function(x) as.numeric(x) < 70)
##################
 
p <- flip(p, 73, 85)
p <- p +
  geom_facet(panel = "Ocean Basin", 
             data = Broad_strata, 
             geom = geom_col, 
             #geom = geom_col_pattern, # USE THIS IF SOME STRATA ARE CROSS-HATCHED
             mapping = aes(x = proportion, 
                           #pattern = pattern_type,
                           fill = Broad), 
             position = 'stack',
             orientation = 'y',
             width = 1) +
  scale_fill_manual(values = broad_stratum_colors) +
#  scale_pattern_manual(values = c('none' = 'none', 'crosshatch' = 'stripe'), guide = 'none') +
  theme(strip.text = element_blank(), 
        element_text(size = 14))   # Remove panel labels

p <- facet_widths(p, c(9, 1))
#jpeg(file = 'results-raw/viridis_tree.jpg', width = 1000, height = 1000)
pdf(file = 'results-raw/viridis_tree_logMut_logMRCA.pdf', width = 7, height = 7)
p
dev.off()

