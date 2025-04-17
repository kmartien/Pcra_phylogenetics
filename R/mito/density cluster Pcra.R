library(strataG)
library(tidyverse)
#library(gdata)
library(densityClust)
library(ape)
#library(ggbiplot)
library(grid)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)

load('data/Pcra.mito.gtype.rda')

dloop.dist <- dist.dna(
  getSequences(mito.g, as.haplotypes = TRUE), 
  model = "TN93"
)
dc.est <- estimateDc(dloop.dist)
dloop.dc <- densityClust(dloop.dist, dc = dc.est)
dloop.dc <- findClusters(dloop.dc)
centers <- names(dloop.dc$rho)[dloop.dc$peaks]
denClust.results <- bind_cols(haplotype = names(dloop.dc$rho), 
               rho = dloop.dc$rho, 
               delta = dloop.dc$delta, 
               halo = dloop.dc$halo, 
               cluster = as.factor(dloop.dc$clusters)) |> 
  mutate(prod.rho.delta = rho * delta,
         cluster.center = ifelse(haplotype %in% centers, TRUE, FALSE)) 
write.csv(denClust.results, file = "clusters.csv", row.names = FALSE)

############################################
###
### Customized version of plotDensityClust

plotDensityClustKKM <-
function (x, type = "all", n = 20, mds = NULL, dim.x = 1, dim.y = 2, 
          col = NULL, alpha = 0.8) 
{
  type <- tolower(type)
  if (any(pmatch(type, "all", nomatch = 0))) 
    type <- c("dg", "gg", "mds")
  df <- data.frame(hap = names(x$rho), rho = x$rho, delta = x$delta, gamma = x$rho * 
                     x$delta, peaks = FALSE, cluster = factor(x$clusters), 
                   halo = x$halo)
  df$peaks[x$peaks] <- TRUE
  if (is.null(col)) {
    num.cols <- max(nlevels(df$cluster), 3)
    col <- if (num.cols <= 8) {
      brewer.pal(num.cols, "Set2")
    }
    else if (num.cols <= 12) {
      brewer.pal(num.cols, "Set3")
    }
    else rainbow(num.cols + 1)[1:num.cols]
  }
  plots <- list(dg = NULL, gg = NULL, mds = NULL)
  if (any(pmatch(type, "dg", nomatch = 0))) {
    plots$dg <- ggplot(df, aes_string(x = "rho", y = "delta"))
    ######################################################################
    # removing the lines that identify the threshold values of rho and delta, since
    # they're never discussed and clutter up the plot
    # if (!any(is.na(x$threshold))) {
    #   rho <- x$threshold["rho"]
    #   delta <- x$threshold["delta"]
    #   thresh.df <- data.frame(x = c(rho, rho), y = c(delta, 
    #                                                  delta), xend = c(rho, Inf), yend = c(Inf, delta))
    #   plots$dg <- plots$dg + geom_segment(aes_string(x = "x", 
    #                                                  xend = "xend", y = "y", yend = "yend"), data = thresh.df, 
    #                                       inherit.aes = F, lineend = "butt")
    # }
    ######################################################################
    if (any(df$peaks)) {
      plots$dg <- plots$dg + geom_label_repel(aes_string(label = "hap",
                                                   color = "cluster"), data = df[df$peaks, ], fontface = "bold",
                                        alpha = alpha) + scale_color_manual(values = col)
      plots$dg <- plots$dg + geom_point(data = df[df$peaks,], size = 3, color = 'gray50', shape = 17)# + scale_color_manual(values = col)
    }
    plots$dg <- plots$dg + geom_point(data = df[!df$peaks,], size = 3, color = "gray50", alpha = alpha) + labs(main = "Decision Graph", x = expression(rho), 
                                                         y = expression(delta), color = "Cluster") + theme(legend.position = "none")
  }
  if (any(pmatch(type, "gg", nomatch = 0))) {
    gg.df <- df[order(df$gamma, decreasing = TRUE), ]
    gg.df <- gg.df[1:n, , drop = FALSE]
    gg.df$Sample <- 1:nrow(gg.df)
    plots$gg <- ggplot(gg.df, aes_string(x = "Sample", y = "gamma")) + 
      geom_line()
    if (any(gg.df$peaks)) {
      plots$gg <- plots$gg + geom_label(aes_string(label = "hap", 
                                                   color = "cluster"), data = gg.df[gg.df$peaks, 
                                                                                    , drop = FALSE], fontface = "bold", alpha = alpha) + 
        scale_color_manual(values = col)
    }
    plots$gg <- plots$gg + geom_point(data = gg.df[!gg.df$peaks, 
                                                   , drop = FALSE], size = 3, color = "gray50") + labs(y = expression(gamma), 
                                                                                                       color = "Cluster") + theme(legend.position = "none")
  }
  if (any(pmatch(type, "mds", nomatch = 0))) {
    if (is.null(mds)) 
      mds <- cmdscale(x$distance, k = max(dim.x, dim.y))
    
    
    df$x <- mds[, dim.x]
    df$y <- mds[, dim.y]
    plots$mds <- ggplot()
    plots$mds <- if (all(is.na(df$cluster))) {
      plots$mds + geom_point(aes_string(x = "x", y = "y"), 
                             data = df, size = 3, color = "gray50", alpha = alpha)
    }
    else {
      plots$mds + geom_point(aes_string(x = "x", y = "y", 
                                        color = "cluster"), data = df[df$halo, , drop = FALSE], 
                             shape = 21, size = 3) + 
        geom_point(aes_string(x = "x", y = "y", color = "cluster"), 
                   data = df[!df$halo, , drop = FALSE], size = 3, alpha = alpha) + 
        geom_label_repel(aes_string(x = "x", y = "y", label = "hap", color = "cluster"), 
                         data = df[df$peaks, , drop = FALSE], size = 6, 
                         fontface = "bold", alpha = alpha) + scale_color_manual(values = col, 
                                                                                na.value = "gray50")
    }
    plots$mds <- plots$mds + labs(x = paste("Dimension", 
                                            dim.x), y = paste("Dimension", dim.y)) + theme(legend.position = "none")
  }
  has.plot <- !sapply(plots, is.null)
  switch(sum(has.plot), print(plots[[which(has.plot)]]), {
    plots <- plots[has.plot]
    if ("mds" %in% names(plots)) plots$nrow <- 2 else plots$ncol <- 2
    do.call(grid.arrange, plots)
  }, {
    plots$layout_matrix <- matrix(c(1, 3, 2, 3), nrow = 2)
    do.call(grid.arrange, plots)
  })
}
### End plotDensityClustKKM
############################################

png(filename = 'densityClust.plot.png', width = 1500, height = 1500, res = 200)
plotDensityClustKKM(dloop.dc)
dev.off()
