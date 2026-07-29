library(strataG)
library(swfscMisc)
library(randomForest)

scheme <- "NAtl_vs_Other"

###############################################
# Mitogenomes

load('data/Pcra.mito.gtype.rda')

g <- stratify(mito.g, scheme)
g@description <- paste0('Pcra.mito.', scheme)

stratified.rfDF <- gtypes2rfDF(g)

# Identify groups of sites that are perfectly correlated and remove all but one from each group
# stratified.rfDF.num <- do.call("cbind",lapply(stratified.rfDF[-1],function(x){as.numeric(x)}))
# site.correl <- cor(stratified.rfDF.num)
# 
# to.remove <- sort(unique(do.call('c',lapply(1:dim(site.correl)[1],function(x){
#   matches <- which(abs(site.correl[x,])>=0.99)
#   if (length(matches)>1) return(matches[1:(length(matches)-1)]) else return()
# }))))
# stratified.rfDF <- stratified.rfDF[,-(to.remove+1)]

#making sure correlated sites were removed correctly
#sites.kept <- sapply(colnames(stratified.rfDF)[-1],function(x){as.numeric(substr(x,start=6,stop=nchar(x)))})
#seq.mat <- as.matrix(seqs)
#write.fasta(seq.mat[,sites.kept],file=paste(description,"seqs.fast",sep=""))

freq <- table(stratified.rfDF$stratum)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

rf <- randomForest(
  stratum ~ .,
  stratified.rfDF,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000
)

###############################################
# Mitogenomes

load('data/Pcra.CR.gtype.rda')

g <- stratify(CR.g, scheme)
g@description <- paste0('Pcra.CR.', scheme)

stratified.rfDF <- gtypes2rfDF(g)

# Identify groups of sites that are perfectly correlated and remove all but one from each group
# stratified.rfDF.num <- do.call("cbind",lapply(stratified.rfDF[-1],function(x){as.numeric(x)}))
# site.correl <- cor(stratified.rfDF.num)
# 
# to.remove <- sort(unique(do.call('c',lapply(1:dim(site.correl)[1],function(x){
#   matches <- which(abs(site.correl[x,])>=0.99)
#   if (length(matches)>1) return(matches[1:(length(matches)-1)]) else return()
# }))))
# stratified.rfDF <- stratified.rfDF[,-(to.remove+1)]

#making sure correlated sites were removed correctly
#sites.kept <- sapply(colnames(stratified.rfDF)[-1],function(x){as.numeric(substr(x,start=6,stop=nchar(x)))})
#seq.mat <- as.matrix(seqs)
#write.fasta(seq.mat[,sites.kept],file=paste(description,"seqs.fast",sep=""))

freq <- table(stratified.rfDF$stratum)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

rf <- randomForest(
  stratum ~ .,
  stratified.rfDF,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000
)

