#' sere
#'
#' SERE [Single-parameter quality control and sample comparison for RNA-Seq]
#'
#' Code based on https://github.com/eigenv/SERE
#'
#' obs.count: table of gene and sample level read counts
#'
#' min.reads: minimum number of reads for a gene for exclusion
#'
#' @return
#' @export
#'
#' @examples
#'
sere <- function(obs.count, min.reads = 1) {

  num.samp <- ncol(obs.count)

  tot.count <- sum(obs.count)
  samp.count <- colSums(obs.count)

  row.count <- rowSums(obs.count)

  idx <- which(row.count > min.reads)
  obs.count <- obs.count[idx, ]
  row.count <- row.count[idx]
  num.genes <- nrow(obs.count)

  expt.count <- matrix(NA, nrow = num.genes, ncol = num.samp)

  for(gene.idx in seq(num.genes)) {
    expt.count[gene.idx, ] <- samp.count * row.count[gene.idx] / tot.count
  }

  disp.sum <- sum((obs.count - expt.count) ^ 2 / expt.count)
  sere <- sqrt(disp.sum / (num.genes * (num.samp - 1)))

  return(sere)
}