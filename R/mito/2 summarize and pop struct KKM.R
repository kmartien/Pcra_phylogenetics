library(strataG)
library(hfufs)
library(tidyverse)
source("R/functions/summarizeLoci.bystrata.R")  ### this can go away once Eric fixes the by.strata bug

load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('../Pcra.database.data/data/ASNGS.rda')

scheme <- "mito.Broad"
n.reps.pvals <- 10000  ### change this to a small number for debugging purposes

##############################################################################
## Mitogenome diversity and differentiation

load('data/Pcra.mito.gtype.rda')

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
names(strata.smry)[9:10] <- c("Nei_pi_nuc_div","other_mean.nuc.diversity")

# overall summaries
overall.smry <- summarizeLoci(g)
overall.smry$mean.within.nuc.divergence <- NA
overall.smry$meannuc.diversity <- mean(nucleotideDiversity(g), na.rm = TRUE)
names(overall.smry) <- names(strata.smry)
overall.smry$stratum <- "Overall"

mito.smry.tbl <- rbind(strata.smry, overall.smry)
mito.smry.tbl[6:10] <- lapply(mito.smry.tbl[6:10], signif, digits = 3)

hap.freq.strata <- as.matrix(alleleFreqs(g,by.strata = TRUE)$mito.hap)
hap.freq <- as.data.frame.matrix(cbind(hap.freq.strata, alleleFreqs(g)$mito.hap)) |> 
  rownames_to_column(var = 'mito.hap')

mito_to_CR_key <- ASNGS |> 
  left_join(mito.haps) |> 
  left_join(CR.haps) |> 
  select(c(mito.hap, CR.hap)) |> 
  distinct()
hap.freq <- right_join(mito_to_CR_key, hap.freq) |> 
  mutate(CR.hap = paste0("Hap_", CR.hap))

# pairwise comparisons
# pairwiseTest fails if you send it strata with only one sample, so getting rid of those for the differentiation/divergence tests
g.5 <- g[,,which(getNumInd(g, by.strata = TRUE)$num.ind > 5)]
pws <- pairwiseTest(g.5, nrep = n.reps.pvals, model = mdl)
pws.struct <- pairwiseSummary(pws)
nuc.dvg <- nucleotideDivergence(g.5, model=mdl)

pws.results <- merge(
  pws.struct,
  nuc.dvg$between[, c("strata.1", "strata.2", "dA", "mean")],
  by = c("strata.1", "strata.2"),
  all.x = TRUE
)
PhiST_mat <- pairwiseMatrix(pws, stat = "PHIst")

overall.struct <- overallTest(g.1, nrep = n.reps.pvals, model = mdl)$result

# output results
label <- getDescription(g)
hap.freq.fname <- paste0('results-raw/', label, "_hapfreq.csv")
smry.fname <- paste0('results-raw/', label, "_diversity.stats.csv")
ovl.fname <- paste0('results-raw/', label, "_overall.pop.struct.csv")
pws.fname <- paste0('results-raw/', label, "_pairwise.pop.struct.csv")
pairwise.mat.fname <- paste0('results-raw/', label, "_PhiST_mat.csv")

write.csv(hap.freq, file = hap.freq.fname)
write.csv(mito.smry.tbl, file = smry.fname)
write.csv(overall.struct, file = ovl.fname)
write.csv(pws.results, file = pws.fname)
write.csv(PhiST_mat, file = pairwise.mat.fname)

# dA for NAtl vs. other
g <- stratify(mito.g, scheme = 'NAtl_vs_Other')
g@description <- paste0('Pcra.mito.', scheme)

taxonomic_status <- pairwiseSummary(pairwiseTest(g, nrep = n.reps.pvals, model = mdl))
nuc.dvg <- nucleotideDivergence(g, model=mdl)

taxonomic_status <- merge(
  taxonomic_status,
  nuc.dvg$between[, c("strata.1", "strata.2", "dA", "mean")],
  by = c("strata.1", "strata.2"),
  all.x = TRUE
)
write.csv(taxonomic_status, 
          file = paste0('results-raw/', label, "_taxonomic_status.csv"))

##############################################################################
## Control region diversity and differentiation

load('data/Pcra.CR.gtype.rda')

# scheme <- "mito.Broad"
# n.reps.pvals <- 10  ### change this to a small number for debugging purposes

g <- stratify(CR.g, scheme)
g@description <- paste0('Pcra.CR.', scheme)

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
names(strata.smry)[9:10] <- c("Nei_pi_nuc_div","other_mean.nuc.diversity")

hap.freq.strata <- alleleFreqs(g,by.strata = TRUE)$CR.hap
hap.freq <- cbind(hap.freq.strata, alleleFreqs(g)$CR.hap)

# overall summaries
overall.smry <- summarizeLoci(g)
overall.smry$mean.within.nuc.divergence <- NA
overall.smry$meannuc.diversity <- mean(nucleotideDiversity(g), na.rm = TRUE)
names(overall.smry) <- names(strata.smry)
overall.smry$stratum <- "Overall"

CR.smry.tbl <- rbind(strata.smry, overall.smry)
CR.smry.tbl[6:10] <- lapply(smry.tbl[6:10], signif, digits = 3)
#write.csv(smry.tbl, file = "results-raw/diversity_summary_table.csv")

# pairwise comparisons
# pairwiseTest fails if you send it strata with only one sample, so getting rid of those for the differentiation/divergence tests
g.5 <- g[,,which(getNumInd(g, by.strata = TRUE)$num.ind > 5)]
pws <- pairwiseTest(g.5, nrep = n.reps.pvals, model = mdl)
pws.struct <- pairwiseSummary(pws)
nuc.dvg <- nucleotideDivergence(g.5, model=mdl)

pws.results <- merge(
  pws.struct,
  nuc.dvg$between[, c("strata.1", "strata.2", "dA", "mean")],
  by = c("strata.1", "strata.2"),
  all.x = TRUE
)
PhiST_mat <- pairwiseMatrix(pws, stat = "PHIst")

overall.struct <- overallTest(g.5, nrep = n.reps.pvals, model = mdl)$result

# output results
label <- getDescription(g)
hap.freq.fname <- paste0('results-raw/', label, "_hapfreq.csv")
smry.fname <- paste0('results-raw/', label, "_diversity.stats.csv")
ovl.fname <- paste0('results-raw/', label, "_overall.pop.struct.csv")
pws.fname <- paste0('results-raw/', label, "_pairwise.pop.struct.csv")
pairwise.mat.fname <- paste0('results-raw/', label, "_PhiST_mat.csv")

write.csv(hap.freq, file = hap.freq.fname)
write.csv(CR.smry.tbl, file = smry.fname)
write.csv(overall.struct, file = ovl.fname)
write.csv(pws.results, file = pws.fname)
write.csv(PhiST_mat, file = pairwise.mat.fname)

######################################################################
# dA for NAtl vs. other
g <- stratify(CR.g, scheme = 'NAtl_vs_Other')
g@description <- paste0('Pcra.CR.', scheme)

taxonomic_status <- pairwiseSummary(pairwiseTest(g, nrep = n.reps.pvals, model = mdl))
nuc.dvg <- nucleotideDivergence(g, model=mdl)

taxonomic_status <- merge(
  taxonomic_status,
  nuc.dvg$between[, c("strata.1", "strata.2", "dA", "mean")],
  by = c("strata.1", "strata.2"),
  all.x = TRUE
)
write.csv(taxonomic_status, 
          file = paste0('results-raw/', label, "_taxonomic_status.csv"))

########################################################################
# Tajima's D and Fus Fs

CR_Tajima_D <- tajimasD(CR.g)
mito_Tajima_D <- tajimasD(mito.g)

CR_overall_smry <- summarizeLoci(CR.g)
CR_Fu <- afufs(n = CR_overall_smry$num.genotyped, k = CR_overall_smry$num.haplotypes, theta = CR_overall_smry$haplotypic.diversity)
mito_overall_smry <- summarizeLoci(mito.g)
mito_Fu <- afufs(n = mito_overall_smry$num.genotyped, k = mito_overall_smry$num.haplotypes, theta = mito_overall_smry$haplotypic.diversity)
