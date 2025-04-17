library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

tree2 <- read.beast(file = 'BEAST/xml/aln.3.CDS.rRNA.best-MCC.tree')
clades <- tibble(
  clade = c(1,3:8),
  tip1 = c('Pcra.mito.55','Pcra.mito.22','Pcra.mito.09','Pcra.mito.57','Pcra.mito.67','Pcra.mito.66','Pcra.mito.39'),
  tip2 = c('Pcra.mito.10','Pcra.mito.03','Pcra.mito.45','Pcra.mito.33','Pcra.mito.20','Pcra.mito.07','Pcra.mito.01')
)
clade_members <- lapply(1:nrow(clades), function(x){
  node <- MRCA(as_tibble(tree2), clades$tip1[x], clades$tip2[x])$node
  offspring(as_tibble(tree2), node) |> 
    filter(!is.na(label)) |> 
    select(label) |> 
    mutate(clade = clades$clade[x]) |> 
    rename(mito.hap = label)
}) |> bind_rows() |> 
  bind_rows(tibble(mito.hap = 'Pcra.mito.28', clade = 2))

save(clade_members, file = 'data/clade_membership.rda')
