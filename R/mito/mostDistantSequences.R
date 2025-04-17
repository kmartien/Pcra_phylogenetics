library(tidyverse)
library(apex)

normie.haps <- apex::read.multiFASTA("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/Pcra.normie.mitos.fasta")
unique.haps <- apex::read.multiFASTA("alignments/Pcra.unique.mitos.lt6.Ns.fasta")
num.seqs <- 2

chosen.seqs <- mostDistantSequences(unique.haps, num.seqs)
write.csv(chosen.seqs, file = "results-raw/most.distant.haps.csv")
#' @title Most Distant Sequences
#' @description Finds the set of sequences that are on the edges of the cloud of
#'   distances. These are the ones that have the greatest mean pairwise distance
#'   and greatest variance in distances.
#'
#' @param x a set of sequences or a \linkS4class{gtypes} object with sequences.
#' @param num.seqs number of sequences to return. If \code{NULL} (default), all
#'   sequences are returned from most to least distant.
#' @param model a character string specifying the evolutionary model to be used.
#'   See \link{dist.dna} for more information.
#' @param pairwise.deletion a logical indicating whether to delete sites with
#'   missing data. See \link{dist.dna} for more information.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#'
#' @return a vector of the \code{num.seqs} sequence names that are the most
#'   divergent sorted from greatest to least distant.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.haps)
#' 
#' mostDistantSequences(dolph.haps, 5)
#' 
#' @export
#' 
mostDistantSequences <- function(
  x, num.seqs = NULL, model = "raw", pairwise.deletion = TRUE, simplify = TRUE
) { 
  
  result <- sapply(
    apex::getSequences(x, simplify = FALSE),
    function(dna) {
      # calculate distance between sequences
      seq.dist <- ape::dist.dna(
        dna, 
        model = model, 
        pairwise.deletion = pairwise.deletion, 
        as.matrix = TRUE
      )
      
      # initialize with furthest sequence
      opt <- options(warn = -1)
      seq.cmd <- stats::cmdscale(seq.dist, k = length(dna) - 1)
      options(opt)
      range.mean <- apply(seq.cmd, 2, function(p) {
        p %>% 
          range(na.rm = TRUE) %>% 
          mean()
      })
      dist.to.centroid <- (t(seq.cmd) - range.mean) %>% 
        t() %>% 
        apply(1, function(p) {
          sqrt(sum(p ^ 2, na.rm = TRUE))
        }) %>% 
        sort(decreasing = TRUE)
      
      num.seqs <- if(is.null(num.seqs)) {
        length(dna) 
      } else {
        min(num.seqs, length(dna))
      }
      ids <- rep(as.character(NA), length = num.seqs)
      ids[1] <- names(dist.to.centroid)[1]
      if(num.seqs == 1) return(ids)
      
      for(i in 2:length(ids)) {      
        current.ids <- ids[1:(i - 1)]
        # add sequence with greatest mean distance and variance to current set
        ids[i] <- setdiff(rownames(seq.dist), current.ids) %>% 
          purrr::map(function(id) {
            dist.to.current <- seq.dist[id, current.ids]
            tibble::tibble(
              id = id,
              mean = mean(dist.to.current, na.rm = TRUE),
              var = stats::var(dist.to.current, na.rm = TRUE)
            )
          }) %>% 
          dplyr::bind_rows() %>% 
          dplyr::arrange(dplyr::desc(.data$mean), dplyr::desc(.data$var)) %>% 
          dplyr::slice(1) %>% 
          dplyr::pull(.data$id)
      }
      ids
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}