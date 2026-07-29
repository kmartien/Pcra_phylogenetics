# plot PSMC plots with bootstraps, using output from psmc_plot.pl
# modified by K. Hernandez May 2025
# Updated Feb 2026 for visual consistency across scripts

library(viridis)

# Disable scientific notation globally for this R session
options(scipen = 999)

PSMC_folder <- "hPSMC_Dsuite/hPSMC"
allfiles = list.files(path = PSMC_folder, pattern="txt")

pop1 <- "Pcra_z0018462"
pop2 <- "Pcra_z0027510"
pop3 <- "Pcra_z0045928"
pop4 <- "hPSMC_z0018462_z0027510"
pop5 <- "hPSMC_z0018462_z0045928"
pop6 <- "hPSMC_z0027510_z0045928"

leg1 <- c("ETP", "ATL", "MHI", "hPSMC:ETP-ATL", "hPSMC:ETP-MHI", "hPSMC:ATL-MHI") 

# Save as PDF
pdf(paste0("results-raw/Pseudorca", "_psmc_plot.pdf"), width = 6, height = 3.5)

### Set min and max values for plot axes (True Absolute Ne Scale)
xmin = 0
xmax = 225000
ymin = 0
ymax = 180000 # Absolute scale ceiling

### Set line colors 
c1 = c(229, 198, 27)/255
c2 = c(72, 21, 104)/255
c3 = c(117, 117, 117)/255
c4 = c(0, 0, 0)/255
c5 = c(0, 0, 0)/255
c6 = c(0, 0, 0)/255

mycols1 = c(rgb(c1[1], c1[2], c1[3], alpha=1), rgb(c2[1], c2[2], c2[3], alpha=1),
            rgb(c3[1], c3[2], c3[3], alpha=1), rgb(c4[1], c4[2], c4[3], alpha=1),
            rgb(c5[1], c5[2], c5[3], alpha=1), rgb(c6[1], c6[2], c6[3], alpha=1))

transp = 0.05 
mycols2 = c(rgb(c1[1], c1[2], c1[3], alpha=transp), rgb(c2[1], c2[2], c2[3], alpha=transp),
            rgb(c3[1], c3[2], c3[3], alpha=transp), rgb(c4[1], c4[2], c4[3], alpha=0),
            rgb(c5[1], c5[2], c5[3], alpha=0),     rgb(c6[1], c6[2], c6[3], alpha=0))

### Generate canvas with clean styling configuration
# mar[2]=5.5 clears space for titles; cex.axis=0.8 and cex.lab=1.0 make titles slightly larger than ticks
par(mar=c(4, 6.8, 1, 1), family = "sans", cex = 0.8, cex.axis = 0.8, cex.lab = 1.0) 

# Initialize empty plot without automatic boxes or scales
plot(1, 1, type="n", axes=FALSE, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")

# INJECT BACKGROUND GRID: Wider vertical light-gray lines (lwd = 1)
abline(v = seq(0, 200000, by = 50000), col = "gray90", lty = 1, lwd = 1)

# Axis titles (line=4.0 pushes the y-axis title out past the large numbers)
title(xlab="Years before present", line=2.5)
title(ylab=expression("Nuclear Effective Population Size (" * italic(N)[e] * ")"), line=4.0)

# ==============================================================================
# VISUAL CONSISTENCY REVISIONS (AXIS MAPPING)
# ==============================================================================
# 1. Thicker standard continuous black L-shaped axis lines
box(bty = "l", lwd = 1.2, col = "black")

# 2. X-axis: labels every 50,000 (no tick marks)
at.x <- seq(0, 200000, by = 50000)
lab.x <- prettyNum(at.x, big.mark = ",")
axis(1, at = at.x, labels = lab.x, tick = FALSE, las = 1)

# 3. Y-axis: micro-ticks every 10,000 (matching tick thickness to thicker line)
axis(2, at = seq(0, ymax, by = 10000), labels = FALSE, lwd = 0, lwd.ticks = 1.0)

# 4. Y-axis: macro-labels every 50,000
at.y_labels <- seq(0, ymax, by = 50000)
lab.y <- prettyNum(at.y_labels, big.mark = ",")
axis(2, at = at.y_labels, labels = lab.y, tick = FALSE, las = 1)
# ==============================================================================

legend("topright", lwd=3, col=mycols1, legend=leg1, bty="n", lty = c(1, 1, 1, 1, 2, 3))

### Function to add plot lines for each sample (Modified to process absolute Ne)
psmc_plot_fill = function(){
  dfiles = paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
  for (i in 2:length(dfiles)){
    bb = read.table(dfiles[i])
    bb$V1[aa$V1 < 10000] <- 10000 # truncating at 10,000 yr bp as is standard for PSMC
    bb$V2 <- bb$V2 * 10000 # Convert to absolute scale
    lines(bb$V1, bb$V2, type="s", col=mycols2[nn], lwd=1)
  }
  aa = read.table(dfiles[1])
  aa$V1[aa$V1 < 10000] <- 10000 # truncating at 10,000 yr bp as is standard for PSMC
  aa$V2 <- aa$V2 * 10000 # Convert to absolute scale
  lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2)
}

# Plot Main Lineages & Bootstraps
samplename=pop1; nn=1; psmc_plot_fill()
samplename=pop2; nn=2; psmc_plot_fill()
samplename=pop3; nn=3; psmc_plot_fill()

# Plot hPSMC Data Spans
samplename=pop4; nn=4
dfiles = paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa = read.table(dfiles[1]); aa$V2 <- aa$V2 * 10000
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 1)

samplename=pop5; nn=5
dfiles = paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa = read.table(dfiles[1]); aa$V2 <- aa$V2 * 10000
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 2)

samplename=pop6; nn=6
dfiles = paste0(PSMC_folder, "/", allfiles[grep(pattern=samplename, allfiles)])
aa = read.table(dfiles[1]); aa$V2 <- aa$V2 * 10000
lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2, lty = 3)

dev.off()