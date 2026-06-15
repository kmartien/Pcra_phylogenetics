library(tidyverse)

plup <- read.csv("data-raw/Pcra-Allruns-edit.dist4.pileup.ref.csv")

pool_prob <- plup %>%
  mutate(
    reads_A = A * coverage,
    reads_C = C * coverage,
    reads_G = G * coverage,
    reads_T = T * coverage,
    reads_N = N * coverage,
    reads_X = X. * coverage
  ) |> 
  group_by(ref.pos) %>%
  summarise(
    total_coverage = sum(coverage),
    freq_A = sum(reads_A)/total_coverage,
    freq_C = sum(reads_C)/total_coverage,
    freq_G = sum(reads_G)/total_coverage,
    freq_T = sum(reads_T)/total_coverage,
    freq_N = sum(reads_N)/total_coverage,
    freq_X = sum(reads_X)/total_coverage,
    n_individuals = n(),
    .groups = "drop"
  )

rare_alleles <- plup |> 
  filter(
    (A > 0.2 & A < 0.5) |
      (C > 0.2 & C < 0.5) |
      (G > 0.2 & G < 0.5) |
      (T > 0.2 & T < 0.5)
  )

Falklands_SNP <- plup |> 
  filter(
    (id %in% c('z0114824', 'z0114825', 'z0114827', 'z0114829', 'z0114830', 'z0114836', 'z0114838', 'z0114840', 'z0114842')) &
      (ref.pos %in% 16250:16275)
  )
# Falklands SNP is at pos 16269

pool_prob_Falklands_SNP <- filter(pool_prob, ref.pos == 16269)
Falklands_SNP <- plup |> 
  filter(ref.pos == 16269)
