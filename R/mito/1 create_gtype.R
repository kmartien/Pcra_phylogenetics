library(ape)
library(strataG)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))
unique.mito.seqs <- read.FASTA('alignments/Pcra.unique.mitos.fasta')

mito.strata <- left_join(
  select(mito.haps, Animal.ID),
  CR.haps,) |> 
  left_join(Pcra.strata) |> 
  rename(id = Animal.ID)

mito.g <- df2gtypes(mito.haps, ploidy = 1, id.col = 1, strata.col = NULL, 
                    loc.col = 2, schemes = mito.strata, 
                    sequences = unique.mito.seqs, description = 'Pcra.mito')

save(mito.g, file = 'data/Pcra.mito.gtype.rda')
