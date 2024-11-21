library(ape)
library(ggplot2)
library(ggtree)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

tree <- read.nexus(file = 'BEAST/xml/aln.3.CDS.rRNA.best-MCC.tree')

mito.strata <- left_join(mito.haps, CR.haps) |> 
  left_join(Pcra.strata)

broad.labels <- lapply(unique(mito.haps$mito.hap), function (h){
  bind_cols(
    mito.hap = h,
    label = filter(mito.strata, mito.hap == h) |> 
      select(Fine) |> 
      unique() |> 
      paste(sep = ', ')
  )
}) |> bind_rows()

broad.tree <- tree
broad.tree$tip.label <- broad.labels$label[order(broad.labels$mito.hap)]

jpeg(file = 'results-raw/og.tree.jpg',width = 1000, height = 1000)
plot(tree, cex = 1)
dev.off()

jpeg(file = 'results-raw/broad.tree.jpg',width = 1000, height = 1000)
plot(broad.tree, cex = 1)
dev.off()
