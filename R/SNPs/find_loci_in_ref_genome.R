library(tidyverse)
library(seqinr)
library(pegas)

ref_chrom_location <- "data-raw/SNPs/Pcra_reference_genome"
loc_seq_location <- "data-raw/SNPs/Pcra_loc_refs"
#sex_linked_locs <- c("Pcra_AMELX_exon6_Ttru_1080bp.fasta")

loc_seq_names <- list.files(loc_seq_location)
#loc_seq_names <- loc_seq_names[-which(loc_seq_names %in% sex_linked_locs)]

loc_deets <- do.call('rbind', lapply(loc_seq_names, function(f){
  loc_name <- read.fasta(paste0(loc_seq_location, "/", f), whole.header = TRUE) |> 
    names() |> 
    strsplit(split = " ")
  return(
    c(f, 
      substring(loc_name[[1]][7],1, nchar(loc_name[[1]][7])-1),
    strsplit(strsplit(loc_name[[1]][1], split = ':')[[1]][2], '-') |> 
      unlist() 
    )
  )
})) |> data.frame()
names(loc_deets) <- c('fname', 'chr', 'start', 'stop')
loc_deets$start <- as.numeric(loc_deets$start)
loc_deets$stop <- as.numeric(loc_deets$stop)

Pcra_loc_refs <- lapply(unique(loc_deets$chr), function(i){
  ref_chrom <- read.fasta(paste0(ref_chrom_location, '/Pcra_ref_chr', i, '.fasta')) |> 
    unlist() |> 
    unname()
  locs <- filter(loc_deets, chr == i)
  frags <- lapply(1:nrow(locs), function(j){
    loc_seq <- ref_chrom[(locs$start[j]-100):(locs$stop[j] + 100)]
    write(paste0('>Extended_',locs$fname[j]), file = 'data-raw/SNPs/Pcra_extended_loc_refs.fasta', append = TRUE)
    write(paste(toupper(loc_seq), collapse = ""), file = 'data-raw/SNPs/Pcra_extended_loc_refs.fasta', append = TRUE)
  })
  names(frags) <- paste0('Extended_', locs$fname)
})


