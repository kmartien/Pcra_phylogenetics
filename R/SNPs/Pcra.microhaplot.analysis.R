library(tidyverse)
library(vcfR)
library(microhaplot)
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/Pcra.strata.rda")
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data/id.key.rda")

strat <- "Broad"
run.label <- "All_files_aligned_to_extended_Pcra_refs"

# for your dataset: customize the following paths

sam.path <- paste0("data-raw/SNPs/sam.files/", run.label)
label.path <- file.path("data-raw/SNPs/mplot_labels", paste0(run.label, ".label.txt"))
vcf.path <- paste0("vcf/", run.label, ".final.recode.vcf")
out.path <- "results-R/microhaplot"
app.path <- "/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot"

label.dat <- data.frame(fname = list.files(sam.path))
label.dat$LABID <- do.call(c, lapply(strsplit(label.dat$fname, split = "_"), function(x){x[[1]]}))
label.dat <- left_join(
  label.dat,
  (filter(Pcra_id_key, alt.id %in% label.dat$LABID) %>% 
     select(-id.type) %>% 
     left_join(select(Pcra.strata, c(Animal.ID, {{strat}}))) %>% 
     rename(LABID = alt.id)
  )
) %>% 
  select(-Animal.ID) |> 
  mutate(LABID = ifelse(grepl('119L8', fname), paste0(LABID, 'b'), LABID))
write.table(label.dat, file = paste0("data-raw/SNPs/mplot_labels/", run.label,".label.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

vcf <- read.vcfR(vcf.path)
#locus.ignore.file <- read.csv(paste0("microhaplot/",run.label, ".locus_annotation.csv"))

# I've prepped the data, so can just jump straight to running the Shiny app
haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                                  sam.path = sam.path,
                                  out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 8) # use all the cores!

runShinyHaplot(app.path)
