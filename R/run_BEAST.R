library(beautier)
library(beastier)
library(tracerer)
library(ggtree)
library(treeio)
library(ape)
#library(Biostrings)
library(ggplot2)

alignment.label <- "Pcra.mito.unique.aligned2Pele" # or "Pcra.CR.unique"
#Pcra.mito.alignment <- "alignments/Pcra.mito.unique.fasta"
#Pcra.CR.alignment <- "alignments/Pcra.CR.unique.fasta"

models <- list("Pcra.mito.test", "Pcra.CR.test")

input_filename <- paste0(alignment.label, ".fast")
output_filename <- "BEAST/xml/Pcra.mito.test.xml"

create_beast2_input_file(input_filename = input_filename, 
                         output_filename = output_filename,
                         site_model = create_hky_site_model(),
                         clock_model = create_rln_clock_model(),
                         tree_prior = create_yule_tree_prior())

beast2_options <- create_beast2_options(
  input_filename = output_filename
)

# This isn't working - R can't find BEAST
output <- run_beast2_from_options(beast2_options)

# Use tracerer to check effective sample sizes
estimates <- parse_beast_tracelog_file(paste0("BEAST/xml/Pcra.log"))
estimates <- remove_burn_ins(estimates, burn_in_fraction = 0.1)
esses <- calc_esses(estimates, sample_interval = 1000)
table <- t(esses)
colnames(table) <- c("ESS")

# Run TreeAnotator, then...

tree.dat <- read.beast("BEAST/xml/Pcra.mito.final.mcc.tre")
ggtree(tree.dat) + geom_tiplab() + geom_treescale() #+ geom_cladelab(node = c(43, 44), label =  c("60", "75"))
ggtree(tree.dat) + geom_nodelab()
tree <- groupClade(.data = tree.dat, .node=c(155,160))
ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab(aes(subset=(group==0)))

