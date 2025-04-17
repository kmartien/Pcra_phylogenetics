library(tidyverse)
library(ggplot2)
#library(phytools)
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/mito.haps.rda")
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/archive_data.rda")
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/Pcra.strata.rda")
load('data/clade_membership.rda')
#aln.3.tree <- readNexus("BEAST/xml/aln.3-MCC.tree")

#focal.node <- 124

samp.dat <- left_join(mito.haps,
                      (select(archive_data, c(LABID,Latitude, Longitude)) %>% 
                         rename(Animal.ID = LABID))) %>% 
  left_join(Pcra.strata) %>% 
  left_join(clade_members, by = 'mito.hap')
  
# tips <- getDescendants(aln.3.tree, node = focal.node) %>% 
#   as_tibble() %>% 
#   filter(value <= length(aln.3.tree$tip.label))
# tips <- aln.3.tree$tip.label[tips$value]
# 
# #tips <- aln.3.tree$tip.label[which(getDescendants(aln.3.tree, node = focal.node) < length(aln.3.tree$tip.label))]
# samp.dat$clade[which(samp.dat$mito.hap %in% tips)] <- "focal"

jpeg(paste0("results-raw/map_clades.jpg"), height = 1000, width = 1500)
world <- map_data("world")
worldplot <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  #coord_fixed(1.3) + xlim(lon.range) + ylim(lat.range) +
  #  geom_point(data = data_summary, aes(Longitude, Latitude)) +
#  geom_point(data = samp.dat, aes(x=Longitude, y = Latitude, colour = clade, size = 3))
  geom_point(data = filter(samp.dat), aes(x=Longitude, y = Latitude, color = as.factor(clade)), size = 4) 
  # geom_point(data = filter(samp.dat, clade == "focal"), aes(x=Longitude, y = Latitude, color = clade, size = 3)) +
  # scale_colour_manual(name = "clade", values = c("red", "black"))
worldplot
dev.off()
