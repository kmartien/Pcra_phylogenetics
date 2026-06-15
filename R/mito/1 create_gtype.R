library(ape)
library(strataG)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('../Pcra.database.data/data/ASNGS.rda')
load('../Pcra.database.data/data/ASPhyloCR.rda')

## Mitogenome gtypes object #################################################
mito.haps <- left_join(ASNGS, mito.haps)
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

## CR gtypes object ########################################################
CR.haps <- left_join(ASPhyloCR, CR.haps)
unique.CR.seqs <- read.FASTA('alignments/Pcra_unique_CR_haps.fasta')

CR.strata <- left_join(
  select(CR.haps, Animal.ID),
  mito.haps,) |> 
  left_join(Pcra.strata) |> 
  rename(id = Animal.ID)

# add CR.haps for the 7 Aus residents that were sequenced in Australia
CR.haps <- bind_rows(
  CR.haps, 
  tibble(
    Animal.ID = paste0('AusRes', 1:7), 
    CR.hap = 45
  )
)

CR.strata <- bind_rows(
  CR.strata, 
  tibble(
    id = paste0('AusRes', 1:7), 
    mito.hap = NA,
    Location = 'Australia',
    Fine = 'Australia',
    Broad = 'Aus_resident',
    Ocean_Basin = 'Indian Ocean',
    mito.Broad = 'Aus_resident',
    NAtl_vs_Other = 'Other'
  )
)

# Create gtype object
CR.g <- df2gtypes(CR.haps, ploidy = 1, id.col = 1, strata.col = NULL, 
                    loc.col = 2, schemes = CR.strata, 
                    sequences = unique.CR.seqs, description = 'Pcra_CR')

save(CR.g, file = 'data/Pcra.CR.gtype.rda')
