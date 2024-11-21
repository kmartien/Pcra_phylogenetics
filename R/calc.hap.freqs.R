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
  arrange(mito.hap) %>% 
  column_to_rownames("mito.hap")
write.csv(mito.hap.freqs, file = paste0("results-raw/mito.hap.freqs.", strat, ".csv"))

hap.locs <- lapply(1:nrow(mito.hap.freqs), function(i){
  names(mito.hap.freqs)[which(!is.na(mito.hap.freqs[i,]))]
})
names(hap.locs) <- rownames(mito.hap.freqs)
