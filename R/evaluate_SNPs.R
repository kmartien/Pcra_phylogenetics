library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)
library(pegas)
source('R/functions/mplot2tgt.R')

run.label <- "Ten_med_aligned_to_Pcra_ref"
tgt <- mplot2tgt(project = run.label)
vcf <- read.vcfR("vcf/Ten_med_aligned_to_Pcra_ref.SNPs.vcf", convertNA = T)

# THIS IS CURRENTLY A CHIMERA. ABOVE AND BELOW ARE TWO DIFFERENT WAYS TO EVALUATE
# GENOTYPES. I'M RERUNNING THE DE NOVO ASSEMBLY TO GET BETTER ALIGNMENTS, THEN
# WE'LL SEE

# convert vcf to a tibble and add a column 'locus' that combines CHROM and POS
# of the SNPs
tidy.vcf <- vcfR2tidy(vcf, single_frame = TRUE, toss_INFO_column = FALSE,
                      info_fields = c("DP","RO", "AO"), format_fields = c("GT", "RO", "AO"))$dat %>%
  mutate(locus = paste(CHROM, POS, sep= "_")) %>%
  relocate(locus, .after = POS)
