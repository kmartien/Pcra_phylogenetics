library(strataG)
library(tidyverse)
source("R/functions/summarizeLoci.bystrata.R")  ### this can go away once Eric fixes the by.strata bug
load('data/Pcra.mito.gtype.rda')

scheme <- "Ocean_Basin"
n.reps.pvals <- 10  ### change this to a small number for debugging purposes

g <- stratify(mito.g, scheme)
g@description <- paste0('Pcra.mito.', scheme)

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


# by-strata summaries
# strata.smry <- summarizeLoci.bystrata(g) %>% left_join(nucleotideDivergence(g, model=mdl)$within[, c("stratum","mean")])
# ### Once Eric fixes summarizeLoci, replace the line above with this one:
strata.smry <- summarizeLoci(g, by.strata = TRUE) %>% left_join(nucleotideDivergence(g, model=mdl)$within[, c("stratum","mean")])
nuc.dvst.strata <- as.data.frame(sapply(strataSplit(g), function(x) mean(nucleotideDiversity(x), na.rm = TRUE)))
nuc.dvst.strata$stratum <- rownames(nuc.dvst.strata)
strata.smry <- left_join(strata.smry[-1], nuc.dvst.strata)
names(strata.smry)[9:10] <- c("mean.within.nuc.divergence","mean.nuc.diversity")

hap.freq.strata <- alleleFreqs(g,by.strata = TRUE)$mito.hap

# overall summaries
hap.freq <- alleleFreqs(g)$mito.hap
overall.smry <- summarizeLoci(g)
overall.smry$mean.within.nuc.divergence <- NA
overall.smry$meannuc.diversity <- mean(nucleotideDiversity(g), na.rm = TRUE)
names(overall.smry) <- names(strata.smry)
overall.smry$stratum <- "Overall"

smry.tbl <- rbind(strata.smry, overall.smry)

# pairwise comparisons
# pairwiseTest fails if you send it strata with only one sample, so getting rid of those for the differentiation/divergence tests
g.1 <- g[,,which(getNumInd(g, by.strata = TRUE)$num.ind > 1)]
pws.struct <- pairwiseSummary(pairwiseTest(g.1, nrep = n.reps.pvals))
nuc.dvg <- nucleotideDivergence(g.1, model=mdl)

pws.results <- merge(
  pws.struct,
  nuc.dvg$between[, c("strata.1", "strata.2", "dA", "mean")],
  by = c("strata.1", "strata.2"),
  all.x = TRUE
)

overall.struct <- overallTest(g.1, nrep = n.reps.pvals)$result

# output results
label <- paste(getDescription(g), scheme, sep=".")
smry.fname <- paste0('results-raw/', label, ".diversity.stats.csv")
ovl.fname <- paste0('results-raw/', label, ".overall.pop.struct.csv")
pws.fname <- paste0('results-raw/', label, ".pairwise.pop.struct.csv")

write.csv(smry.tbl, file = smry.fname)
write.csv(overall.struct, file = ovl.fname)
write.csv(pws.results, file = pws.fname)
