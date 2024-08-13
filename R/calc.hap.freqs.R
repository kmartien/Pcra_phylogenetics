library(tidyverse)
library(swfscMisc)
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/mito.haps.rda")
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/Pcra.strata.rda")

strat <- "Fine"
mito.samp.dat <- left_join(mito.haps, select(Pcra.strata, c(Animal.ID, {{strat}}))) %>% 
  na.omit() %>% rename(stratum = {{strat}})

mito.hap.freqs <- 
  do.call(bind_rows, lapply(
    unique(mito.samp.dat$stratum), 
    function(i){
      filter(mito.samp.dat, stratum == i) %>% 
        group_by(mito.hap) %>% 
        summarise(freq = n()) %>% 
        bind_cols(stratum = i)
    })) %>% 
  pivot_wider(names_from = stratum, values_from = freq) %>% 
  column_to_rownames("mito.hap")


