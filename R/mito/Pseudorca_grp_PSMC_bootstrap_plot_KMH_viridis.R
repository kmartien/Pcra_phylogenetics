# plot PSMC plots with bootstraps, using output from psmc_plot.pl
#modified by K. Hernandez May 2025

#Update 5 Feb 2026, using viridis colors and similar palette to the mitogenome ms
library(viridis)

PSMC_folder <- "hPSMC_Dsuite/hPSMC"
# step 1 (create the .txt files from .psmc files using psmc_plot.pl) - run in terminal
# step 2 (input .txt files from step 1 to create final plot) - run in R

### Get list of all input files including bootstrap replicates (these are the output text files generated from the utils/psmc_plot.pl -R command)
allfiles=list.files(path = PSMC_folder, pattern="txt")

# get species names based on species name patterns from the filenames
pop1<-"Pcra_z0018462"
pop2<-"Pcra_z0027510"
pop3<-"Pcra_z0045928"
pop4<-"hPSMC_z0018462_z0027510"
pop5<-"hPSMC_z0018462_z0045928"
pop6<-"hPSMC_z0027510_z0045928"
sp7<-""
sp8<-""
sp9<-""

leg1<-c("ETP", "ATL", "MHI", "hPSMC:ETP-ATL", "hPSMC:ETP-MHI", "hPSMC:ATL-MHI") # for legend on plot


# Save as PDF using pdf() and dev.off()
pdf(paste0("results-raw/Pseudorca","_psmc_plot.pdf"),width = 5, height = 3)

### Set min and max values for plot axes
#xmin=1.2e4
#xmax=1e7
xmin <- 0
xmax <- 225000
ymin=0
ymax=18

### Set line colors (first, pick rgb values for each sample, then set main and transparent colors for plot lines)
# color1 (viridis - E5C61BFF - yellow)
c1=c(229, 198, 27)/255

# color 2 (viridis - 481568 - purple)
c2=c(72, 21, 104)/255

# color 3 (gray46 - 757575)
c3=c(117, 117, 117)/255

# color 4 (black - 000000)
c4=c(0, 0, 0)/255

# color5 (black - 000000)
c5=c(0, 0, 0)/255

# color 6 (black - 000000)
c6=c(0, 0, 0)/255



# Main colors (no transparency)
mycols1=c(
  rgb(c1[1], c1[2], c1[3], alpha=1),
  rgb(c2[1], c2[2], c2[3], alpha=1),
  rgb(c3[1], c3[2], c3[3], alpha=1),
  rgb(c4[1], c4[2], c4[3], alpha=1),
  rgb(c5[1], c5[2], c5[3], alpha=1),
  rgb(c6[1], c6[2], c6[3], alpha=1)
)

# Bootstrap replicate colors (with transparency)
transp=0.05 # 0.05; 0 if no bootstraps needed.
mycols2=c(
  rgb(c1[1], c1[2], c1[3], alpha=transp),
  rgb(c2[1], c2[2], c2[3], alpha=transp),
  rgb(c3[1], c3[2], c3[3], alpha=transp),
  rgb(c4[1], c4[2], c4[3], alpha=0),
  rgb(c5[1], c5[2], c5[3], alpha=0),
  rgb(c6[1], c6[2], c6[3], alpha=0)
)

### Generate an empty plot with labeled axes
par(mar=c(3.5,3.75,0.5,0.5))
op <- par(cex = 0.75) # font size

### Generate an empty plot with labeled axes
par(mar=c(3.5,3.75,0.5,0.5))
op <- par(cex = 0.75) # font size

# REMOVED log="x" to make the time axis linear
plot(1, 1, type="n", axes=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")

title(xlab="Years before present", line=2)
title(ylab=expression("Effective population size (x10"^4*")"), line=2.25)

axis(side=2, line=0, labels=F)
axis(side=2, line=-.25, labels=T, tick=F)

# REPLACED the old log-scale loop with clean linear ticks every 50,000 years
at.x <- seq(0, 1000000, by = 50000)
lab.x <- prettyNum(at.x, big.mark = ",")
axis(1, at=at.x, labels=lab.x, las=1)

legend("topright", lwd=3, col=mycols1[c(1,2,3,4,5,6,7,8,9)], legend=leg1, 
       bty="n", lty = c(1, 1, 1, 1, 2, 3, 1, 1, 1))

box()

### Function to add plot lines for each sample
psmc_plot_fill=function(){
  # Get list of input files using "samplename" as the search term
  dfiles=paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
  # Loop through the bootstrap reps (first bootstrap file = the second file from the list above)
  for (i in 2:length(dfiles)){
    bb=read.table(dfiles[i])
    # Plot lines for each bootstrap file using the transparency colors
    lines(bb$V1, bb$V2, type="s", col=mycols2[nn], lwd=1)
  }
  # Read in the first file from the list (this one is for adding the main solid line 
  # on top of the bootstrap lines, given as boot.out.0.txt from the psmc_plot.pl script)
  aa=read.table(dfiles[1])
  # Plot the line using the main color (no transparency)
  lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2)
}

### Plot each sample (nn is the numerical index for the sample - only needs to 
# correspond to the order of samples in the color lists above)
# comment out "psmc_plot_fill()" for unused samplenames.

# samples to plot
# 1
samplename=pop1
nn=1
psmc_plot_fill()
# 2
samplename=pop2
nn=2
psmc_plot_fill()
# 3
samplename=pop3
nn=3
psmc_plot_fill()
# 4. For hPSMC can't use the written function since there are no bootstraps
samplename=pop4
nn=4
dfiles=paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa=read.table(dfiles[1])
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 1)
# 5
samplename=pop5
nn=5
dfiles=paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa=read.table(dfiles[1])
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 2)
# 6
samplename=pop6
nn=6
dfiles=paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa=read.table(dfiles[1])
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 3)

# Close the PDF device
dev.off()

