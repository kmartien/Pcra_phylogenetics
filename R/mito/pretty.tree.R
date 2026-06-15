library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

#tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc_combined/aln.7.CDS-rRNA_constant_orc_combined-MCC.trees')
tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc/aln.3.CDS-rRNA_constant_orc-MCC.tree')

mito.strata <- left_join(mito.haps, CR.haps) |> 
  left_join(Pcra.strata) |> 
  mutate(Broad = ifelse(Broad %in% c('Atl-weird', 'Atl-normal', 'Atl-unknown'), 'N_Atlantic', Broad),
         Broad = ifelse(Broad == 'Atl-South', 'S_Atlantic', Broad),
         Broad = ifelse(Broad == 'Spac', 'S_Pacific', Broad),
         Broad = ifelse(Broad == 'Wpac', 'NPac_Western', Broad),
         Broad = ifelse(Broad == 'Epac', 'NPac_Eastern', Broad),
         Broad = ifelse(Broad == 'CNP', 'NPac_Central', Broad),
         Broad = ifelse(Broad == 'Indian Ocean', 'Indian_Ocean', Broad))

hap.labels <- lapply(unique(mito.haps$mito.hap), function (h){
  bind_cols(
    label = h,
    CR.hap = filter(mito.strata, mito.hap == h) |> 
      select(CR.hap) |> 
      unique() |> 
      paste(sep = ', '),
    fine = filter(mito.strata, mito.hap == h) |> 
      select(Fine) |> 
      unique() |>
      paste(sep = ', '),
    broad = filter(mito.strata, mito.hap == h) |> 
      select(Broad) |> 
      unique() |> 
      paste(sep = ', '),
    OBasin = filter(mito.strata, mito.hap == h) |> 
      select(Ocean_Basin) |> 
      unique() |> 
      paste(sep = ', ')
  )
}) |> bind_rows() |> 
  bind_rows(tibble(label = "ON652887 Pele", CR.hap = 'P. electra')) 

labeled.tree <- full_join(tree, hap.labels, by = 'label') 
labeled.tree@phylo$edge.length[1] <- labeled.tree@phylo$edge.length[1] - labeled.tree@phylo$edge.length[2] + 0.01
labeled.tree@phylo$edge.length[2] <- 0.01
labeled.tree@data$height_0.95_HPD[[68]] <- c(0,0)
x <- groupClade(labeled.tree, MRCA(labeled.tree, 'Pcra.mito.01', 'Pcra.mito.55'))
#x@data$height <- x@data$height*100000
x@phylo$edge.length <- x@phylo$edge.length*1000000
#p <- ggtree(x, aes(linetype = group), mrsd = '2025-01-01') +
p <- ggtree(x, mrsd = '2025-01-01') +
  geom_tree() + theme_tree2() + 
  geom_tiplab(aes(label = fine), size = 3) +
  scale_linetype_manual(values = c(2,1)) +
#  geom_text2(aes(label=node), hjust=1, vjust = -1) +
  geom_point2(aes(subset = posterior > 0.8), size = 2) +
#  geom_text2(aes(label=round((height*1000), 2), subset = posterior > 0.8), vjust = -0.5, hjust = .5, size = 5) +
  scale_x_continuous(breaks = seq(-225000, 0, 50000), minor_breaks = seq(-225000, 0, 10000)) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = 'none') #+
   # geom_text2(aes(label=round((branch), 2)), vjust = -0.5, hjust = .5, size = 5) +
   # geom_text2(aes(label=round(posterior, 2), subset= posterior>0.7), hjust=1, vjust = -1) +
   # geom_range(range = 'height_0.95_HPD', color = 'red', alpha = 0.6, size = 2)
  # geom_hilight(
  # mapping=aes(subset = node %in% c(130, 128, 112, 73, 118, 123, 91, 114),
  #             node = node,
  #             fill = as.factor(node)
  # )
# )
#jpeg(file = 'results-raw/Pcra.mito.tree.CR.hap.jpg', height = 1000, width = 1000)
p
#dev.off()

#nodes2collapse <- c(100, 105)
#collapsed_node_labes <- c('9, 55, 37', '9, 38, 39')
p_collapsed <- p |> collapse(node = 92, mode = 'max', fill = 'coral4') |> #9, 37, 55
  #  collapse(node = 105, mode = 'max', fill = 'darkgreen') |> #9,38,39
  collapse(node = 74, mode = 'max', fill = 'darkred')|> #1, 31
  collapse(node = 82, mode = 'max', fill = 'darkred')|> #45
  collapse(node = 88, mode = 'max', fill = 'darkred')|> #45
  collapse(node = 85, mode = 'max', fill = 'darkred')|> #32
  collapse(node = 81, mode = 'max', fill = 'darkred')|> #2
  collapse(node = 128, mode = 'max', fill = 'purple')|> #26
  #  collapse(node = 93, mode = 'max', fill = 'coral4')|> #30
  collapse(node = 120, mode = 'max', fill = 'darkgreen') |> #50
  collapse(node = 131, mode = 'max', fill = 'coral4') + #23
  geom_text2(aes(subset=(node == 74)), cex=3, vjust=0.2, label="1, 31",hjust = -2) +
  geom_text2(aes(subset=(node == 92)), cex=3, vjust=0.2, label="9, 55",hjust = -1.5) +
  geom_text2(aes(subset=(node == 82)), cex=3, vjust=0.2, label="45",hjust = -1.5) +
  geom_text2(aes(subset=(node == 88)), cex=3, vjust=0.2, label="45",hjust = -3.2) +
  geom_text2(aes(subset=(node == 85)), cex=3, vjust=0.2, label="32",hjust = -1.5) +
  geom_text2(aes(subset=(node == 81)), cex=3, vjust=0.2, label="2",hjust = -1.5) +
  geom_text2(aes(subset=(node == 128)), cex=3, vjust=0.2, label="26",hjust = -1.5) +
  #  geom_text2(aes(subset=(node == 93)), cex=3, vjust=0.2, label="30",hjust = -1.5) +
  geom_text2(aes(subset=(node == 120)), cex=3, vjust=0.2, label="50",hjust = -1.5) +
  geom_text2(aes(subset=(node == 131)), cex=3, vjust=0.2, label="23",hjust = -1.5) 
p_collapsed
#viewClade(p, MRCA(p, 'Pcra.mito.01', 'Pcra.mito.55'))

##############################################################################
# Color-code tips by quasi-ocean basin

Broad_strata <- select(mito.strata, c('mito.hap', 'Broad')) |> 
  distinct() |> 
  rename(label = mito.hap)
p <- ggtree(tree, mrsd = '2025-01-01') +
  geom_tree() + theme_tree2() + 
  #  geom_tiplab(aes(label = fine), size = 3) +
  scale_linetype_manual(values = c(2,1)) +
  #  geom_text2(aes(label=node), hjust=1, vjust = -1) +
  geom_point2(aes(subset = posterior > 0.8), size = 2) +
  #  geom_text2(aes(label=round((height*1000), 2), subset = posterior > 0.8), vjust = -0.5, hjust = .5, size = 5) +
  scale_x_continuous(breaks = seq(-225000, 0, 50000), minor_breaks = seq(-225000, 0, 10000)) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = 'none') #+

p <- p +
  geom_facet(panel = 'Ocean basin', data = Broad_strata, geom = geom_bar, 
             aes(x = label, group = Broad, fill = Broad),
             orientation = 'y', position = 'fill', width = .6) 
p
p <- ggtree(tree) +
  geom_facet(panel = 'Ocean basin', data = Broad_strata, geom = geom_bar, 
             aes(y = mito.hap, group = Broad, fill = Broad),
             position = 'fill', width = .6) +
  scale_fill_manual(values = broad_stratum_colors) 

 
ggplot(Broad_strata) +
  geom_bar(
    aes(y = label, group = Broad, fill = Broad),
    position = 'fill'
  ) +
  scale_fill_manual(values = broad_stratum_colors) 
