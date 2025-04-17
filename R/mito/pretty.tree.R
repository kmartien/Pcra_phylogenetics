library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc/aln.7.CDS-rRNA_constant_orc-MCC.trees')
#tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc/aln.7.CDS-rRNA_constant_orc-Tree.trees')
#blackfish.tree <- read.nexus(file = 'BEAST/xml/aln.6.CDS.rRNA.best_B-combined-MCC.tree')

mito.strata <- left_join(mito.haps, CR.haps) |> 
  left_join(Pcra.strata)


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
p <- ggtree(x, aes(linetype = group)) +
  geom_tree() + theme_tree2() + 
#  geom_tippoint(aes(colour = CR.hap)) + 
  geom_tiplab(aes(label = CR.hap), size = 3) +
  scale_linetype_manual(values = c(2,1)) +
  geom_nodelab(aes(x=branch, label=round((height*1000), 2)), vjust = -0.5, hjust = 1, size = 3) +
  geom_range(range = 'height_0.95_HPD', color = 'red', alpha = 0.6, size = 2)
p
viewClade(p, MRCA(p, 'Pcra.mito.01', 'Pcra.mito.55'))

