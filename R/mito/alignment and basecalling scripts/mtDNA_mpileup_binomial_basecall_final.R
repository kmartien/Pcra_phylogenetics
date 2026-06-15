# Call nucleotide positions in capture-enriched mtDNA based on combined read frequency 
# and binomial probability when "index hopping" causes presence of both allele reads 
# proportional to their frequency in the entire pool.

rm(list = ls())

#source("callBasesPileup_260617_final.R")
source("VanCise_binomial_basecaller_published.R")
require(strataG)

description <- "Pcra-mtdna-NextSeq-PhilParams_20190220_1305" #used for output file names

# read mpileup (.csv) file (output from SWFSC mtDNA pipeline).
plp.fname <- paste(description,"/summaries/","Pcra-mtdna-NextSeq-PhilParams.pileup.ref.csv",sep="")

#define pad length
padlength <- 40
min.cov <- 7

#call bases based on parameter settings
x <- callBasesPileup(plp.fname, min.cov = min.cov, max.prop.only = 7, prop.thresh = 0.85, binom.pr.thresh = 0.95,
                     min.prop = 0.1)
# min.cov = call N if < min.cov
# max.prop.only = minimum coverage needed to call based on frequency
# prop.thresh = Call if frequency of reads >= prop.thresh (and coverage >= max.prop.only)
# binom.pr.thresh = binomial probability threshold (call if >= binom.pr.thresh and frequency >= min.prop)
# min.prop = minimum proportion of all reads per sample for an allele to be considered as potential 
  # allele at a site (needed to exclude very rare reads that are due to genotyping error)

#get sequences
seq.mat <- as.character(x$dna.seqs)

#remove padded ends from sequences
seq.start <- padlength+1
seq.end <- ncol(seq.mat)-padlength
seq.mat<-seq.mat[,c(seq.start:seq.end)]
description <- paste(description,"_nopad")

#write the fasta file
write.fasta(seq.mat, paste(description, "_align.fasta",sep=""))

######SUMMARIES#######
#Count and write file with number of N's per sequence
num.ns <- apply(seq.mat, 1, function(z) sum(tolower(z) == "n"))
num.ns <- data.frame(cbind(id=names(num.ns),total.ns=num.ns))
write.csv(num.ns, paste(description,"_Ns.csv",sep=""))

#write table of variable sites with read counts
varsites <- subset(x$base.freq.wide,variable.calls.at.position==TRUE)
write.csv(varsites, paste(description,"_varsites.csv",sep=""))

#get calls with minor allele
final.call<-x$base.freq.prob.narrow

final.call.sub <- final.call[,c(1,2,4,7:9)]
newtable <- merge(varsites,final.call.sub,by=c('position','id'))
write.csv(newtable,paste(description,"_varsites2.csv",sep=""))

minor.call<-newtable[which(newtable$base==newtable$called.base & newtable$read.prop<0.5),]
write.csv(minor.call,paste(description,"_lowfreq_varsites2.csv", sep=""))

######## THIS DOESN'T WORK CORRECTLY BECAUSE IT IGNORES SITES WITH ZERO COVERAGE ########
# Create a table of samples with number of N's due to low-coverage
#N.lowcov<-newtable[which(newtable$called.base=="N" & newtable$coverage<min.cov & newtable$position>padlength & newtable$position<=seq.end),]
#N.lowcov.unique<-N.lowcov[!duplicated(N.lowcov[,c("position","id")]),]
#low.cov.ns<-data.frame(table(N.lowcov.unique$id))
#colnames(low.cov.ns) <- c("id", "low-cov.Ns")
#write.csv(N.lowcov.unique,paste(description,"_N.locoverage.unique.csv", sep="")) #need to fix to count one N per position/id
#write.csv(N.lowcov,paste(description,"_N.locoverage.csv", sep="")) #need to fix to count one N per position/id
#write.csv(low.cov.ns,paste(description,"_low.cov.n_samples.csv", sep=""))

# create a table of samples with number of N's NOT due to low coverage
N.other<-newtable[which(newtable$called.base=="N" & newtable$coverage>=min.cov & newtable$position>padlength & newtable$position<=seq.end),]
write.csv(N.other,paste(data.set,"/",description,"_N.other.csv", sep=""))
N.other.unique<-N.other[!duplicated(N.other[,c("position","id")]),]
other.ns<-data.frame(table(N.other.unique$id))
colnames(other.ns)<-c("id", "other.Ns")
#write.csv(N.other.unique,paste(description,"_N.other.unique.csv", sep="")) #need to fix to count one N per position/id
#write.csv(other.ns,paste(description,"_other.n_samples.csv", sep=""))

#combined.Ns<-merge(low.cov.ns, other.ns, by = "id", all = TRUE)
combined.Ns<-merge(num.ns, other.ns, by = "id", all = TRUE)
write.csv(combined.Ns, paste(description,"_combined_Ns.csv", sep=""))

#plot frequency of minor allele calls vs. 
# binomial probability of minor allele
hist(minor.call$pr.base)
# frequency of all called bases vs. frequency of reads for major allele 
hist(called.base$read.prop, ylim=c(0,1000))
# need to write pdf file?

#count number of alleles called where the called allele was <50% of the reads (based on probability)
nrow(minor.call)

## Trim 40bp from both ends of fasta file sequences and re-call N's???

