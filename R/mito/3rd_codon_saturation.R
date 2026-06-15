library(strataG)
library(ape)
library(apex)
library(seqinr)
library(tidyverse)
load('data/Pcra.mito.gtype.rda')

I need to redo this so that it only uses the CDS

CDS_BED <- read.csv("data-raw/Pcra_annotated_positions.csv") |> 
  filter(Type=='CDS')
pos2keep <- do.call(c, lapply(1:nrow(CDS_BED), function(i){
  CDS_BED$Start[i]:CDS_BED$Stop[i]
})) |> unique()

g <- mito.g

# choose substitution model
mdl.test.results <- phangorn::modelTest(
  gtypes2phyDat(g), 
  model = c("JC", "F81", "K80", "HKY", "GTR"),
  G = FALSE,
  I = FALSE
) %>% 
  arrange(BIC)

mdl <- switch(
  mdl.test.results$Model[1],
  JC = "JC69",
  F81 = "F81",
  K80 = "K80",
  HKY = "TN93",
  GTR = "TN93"
)

coding_genes <- as.DNAbin(multidna2alignment(g@sequences))[,pos2keep]
distances <- bind_cols(
  corrected = dist.dna(coding_genes, model = 'TN93') |> 
    as.numeric(),
  uncorrected = dist.dna(coding_genes, model = 'raw') |> 
    as.numeric() 
)

regr <- lm(formula = uncorrected ~ corrected, data = distances)

ggplot(distances) +
  geom_point(aes(x=corrected, y=uncorrected))


