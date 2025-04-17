library(dplyr)
library(swfscMisc)

samp.dat <- read.csv("data-raw/ASNGS.haps.csv")
samp.dat$id <- paste0("z0", zero.pad(samp.dat$LabID))

bam.dir <- "data-raw/SNPs/bam.files/aligned.to.CATS.ref"
bam.files <- list.files(path = bam.dir, pattern = ".bam")
bam.files <- do.call(rbind, lapply(bam.files, function(x){
  f.size <- file.info(paste(bam.dir,x, sep = "/")) %>% select(size) 
  id <- substr(x, start = 1, stop = 8)
  return(c(f.name = x, id = id, f.size = f.size$size))
})) %>% data.frame() %>% left_join(samp.dat, by = "id")
bam.files$f.size <- as.numeric(bam.files$f.size)
bam.files <- filter(bam.files, f.size > 7000) #filters out the .bai files

write.csv(bam.files, file = "data-raw/bam.files.size.csv")
