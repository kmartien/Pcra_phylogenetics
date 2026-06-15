library(ape)
library(phytools)
library(treeio)
library(tidyverse)
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/mito.haps.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

mito.haps <- filter(mito.haps, !is.na(mito.hap))

#tree <- read.beast(file = 'BEAST/xml/aln.7.CDS-rRNA_constant_orc_combined/aln.7.CDS-rRNA_constant_orc_combined-MCC.trees')
#tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc/aln.3.CDS-rRNA_constant_orc-MCC.tree')
tree <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc_logMRCA/aln.3.CDS-rRNA_constant_orc_logMRCA-MCC.tree')
tree <- tree@phylo

mito.strata <- left_join(mito.haps, CR.haps) |> 
  left_join(Pcra.strata) 

ob_df <- mito.strata |> 
  select(c(mito.hap, Ocean_Basin)) |> 
  distinct() |> 
  # assign mito hap 18 to either South or North pacific
  filter(!(mito.hap == 'Pcra.mito.18' & Ocean_Basin == 'North Pacific')) |> 
  # remove heteroplasmic haplotypes; they're collapsed to their base hap
  # in the tree
  filter(!mito.hap %in% c('Pcra.mito.68', 'Pcra.mito.69', 'Pcra.mito.70'))

ob_vec <- ob_df$Ocean_Basin
names(ob_vec) <- ob_df$mito.hap

##############################################################################
# Using ace
##############################################################################
pdf("results-raw/ape_ancestral_state_logMRCA.pdf")
plotTree(tree,type="fan",fsize=0.8,ftype="i")

#cols<-setNames(palette()[1:length(unique(ob_vec))],sort(unique(ob_vec)))
#tiplabels(pie=to.matrix(ob_vec,sort(unique(ob_vec))),piecol=cols,cex=0.3)
# add.simmap.legend(colors=cols,
#                   ob_vec=0.9*par()$usr[1],
#                   #y=-max(nodeHeights(tree)),
#                   fsize=0.8)

fitER<-ace(ob_vec,tree,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)

nodelabels(node=1:tree$Nnode+Ntip(tree),
#nodelabels(node=1:tree$Nnode,
                      pie=fitER$lik.anc,piecol=cols,cex=0.5)
#tiplabels(pie=to.matrix(ob_vec,sort(unique(ob_vec))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,
                  prompt=FALSE,
                  ob_vec=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),
                 # x = 0.8,
                  fsize=0.8)
dev.off()

##############################################################################
# Using ancr
##############################################################################

# create a transition matrix defining the state changes that are allowed (i.e.,
# an animal can move from the South Atlantic to the Indian Ocean, but not from
# South Atlantic straight to South Pacific). Use different integers for allowed 
# transitions so that the rates can be estimated independently

ordered_model <- matrix(c(
  0,1,0,0,0,
  2,0,3,4,0,
  0,5,0,6,7,
  0,8,9,0,10,
  0,0,11,12,0),5,5,
  byrow = TRUE,
  dimnames = list(c("North Atlantic", "South Atlantic", "Indian Ocean", "South Pacific", "North Pacific"),
                  c("North Atlantic", "South Atlantic", "Indian Ocean", "South Pacific", "North Pacific")))
ordered_model  

library(foreach)
library(doParallel)
niter<-10 ## set iterations
## set ncores and open cluster
ncores<-min(niter,detectCores()-1)
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
all_fits<-foreach(i=1:niter)%dopar%{ 
  obj<-list()
  class(obj)<-"try-error"
  while(inherits(obj,"try-error")){
    obj<-try(phytools::fitMk(squamate.tree,
                             squamate.toes,model=ordered_model,pi="fitzjohn",
                             logscale=sample(c(FALSE,TRUE),1),
                             opt.method=sample(c("nlminb","optim"),1),
                             rand.start=TRUE))
  }
  obj
}
stopCluster(mc) ## stop cluster
lnL<-sapply(all_fits,logLik)
lnL

