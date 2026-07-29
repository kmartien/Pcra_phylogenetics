library(grid)
library(ggtext)
library(readr)
library(svglite)
library(dplyr)
library(tidyverse)

# Ne used in this simulation was 29,700 (but for the plot scaling use 2.97; it's divided by 10,000) 
# e.g., z0076728.z41749 #2.97 

#Read in the data
species1<-"z0018462"
species2<-"z0045928"
Ne<-1.8988

simulations.combined <- read_table(paste0(species1,"_",species2,"_mut4.90E-10_simulations.combined.txt"), 
                                   col_names = FALSE, col_types = cols(X3 = col_character()))

#Plot all of the simulated data
##!Make sure you change the 2.97 to your Ne. This is important in plotting the 
# shaded areas and drawing the red horizontal lines (indicating the 1.5*Ne and 
# 10*Ne range where you look at how the empirical data is best encased by the simulated data).
dataset <- simulations.combined #250000-600000
dataset$X1[dataset$X1 == 0] <- 1 # get rid of '0' for logarithmic conversion
y1_lim=Ne*1.5
y2_lim=Ne*10

# get time blocks from dataset
blocks<-unique(dataset$X3)
hpsmcNo<- print(length(blocks))
hpsmc<-blocks[hpsmcNo]
blocks<-as.integer(blocks)
blocks<-sort(blocks,decreasing = TRUE)

### plot.margin = margin(top, right, bottom, left)
# change line size for upper/lower CI's (from Cahill et al.: "Divergence is inferred 
# to have occurred between the simulated divergence times of 300–400 ka (red shaded region),
# as these are the closest simulations with transition times that do not intersect the transition
# time of real data.

sp1_sp2 <- ggplot() +
  geom_step(data=dataset[dataset$X3==hpsmc,], aes(x=log10(X1), y=X2, alpha=1, size=0.1)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-30],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-29],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-28],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-27],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-26],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-25],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-24],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-23],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-22],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-21],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-20],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-19],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-18],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-17],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-16],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-15],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-14],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-13],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-12],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  # geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-11],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-10],], aes(x=log10(X1), y=X2, alpha=0.9, size=1.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-9],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-8],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-7],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-6],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-5],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-4],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-3],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-2],], aes(x=log10(X1), y=X2, alpha=0.9, size=0.05)) +
  geom_step(data=dataset[dataset$X3==blocks[hpsmcNo-1],], aes(x=log10(X1), y=X2, alpha=0.9, size=1.05)) +
  theme_classic() +
  labs(
    x = "log(Years before present)",
    y = "Ancestral <b><i>N<sub>e</sub></i></b> (x10^4)"
  ) +
  theme(
    axis.title.y = element_markdown(size = 12, face = "bold") # This ensures the title is bold and 12pt
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  theme(axis.title.x = element_text(face = "bold", size=12)) +
  theme(axis.text.x = element_text(face = "bold", size=12, colour="black")) +
  theme(axis.text.y = element_text(face = "bold", size=12, colour="black")) +
  #theme(axis.title.y = element_text(face = "bold", size=12)) + 
  theme(plot.margin = unit(c(0.0, 0.5, 0.0, 0.0), "in")) +
  coord_cartesian(xlim=c(4, 7), ylim=c(0,50)) + 
  scale_alpha_continuous(range=c(0.3,1), guide=FALSE) +
  scale_size(range=c(0.2,1.5), breaks = c(0.2,1.5), guide = FALSE) +
  labs(colour="Species") +
  geom_hline(yintercept=y1_lim, linetype="dashed", color = "red", linewidth=0.8) +
  geom_hline(yintercept=y2_lim, linetype="dashed", color = "red", linewidth=0.8) +
  # Shade area under y_lim
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = y1_lim),
            alpha = 1/5,
            fill = "grey") +
  # Shade area above y_lim
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = y2_lim, ymax = Inf),
            alpha = 1/5,
            fill = "grey")

pdf(paste0(species1,"_",species2,"_sim_plot.pdf"), width=6, height=3, pointsize=12)
sp1_sp2
dev.off()


