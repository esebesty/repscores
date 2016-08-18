#' npidr
#'
#' npIDR [non-parametric Irreproducible Discovery Rate]
#'
#' Code based on https://github.com/pervouchine/npIDR
#'
#' obs.count: table of gene and sample level read counts
#'
#' pool.type: "sum" or "max"
#'
#' @return
#' @export
#'
#' @examples
#'
npidr <- function(obs.count, bin.size = 1, pool.type = "sum") {

  data <- round(obs.count / bin.size, digits = 0)

  acount1 <- as.data.frame(table(data[, 1]))
  acount2 <- as.data.frame(table(data[, 2]))

  absolute <- merge(acount1, acount2, by = 1, all = T)

  absolute[is.na(absolute)] <- 0

  absolute$sum <- absolute$Freq.x + absolute$Freq.y

  ccount1 <- as.data.frame(table(data[data[, 2] == 0, 1]))
  ccount2 <- as.data.frame(table(data[data[, 1] == 0, 2]))

  conditional <- merge(ccount1, ccount2, by = 1, all = T)

  conditional[is.na(conditional)] <- 0

  conditional$sum = conditional$Freq.x + conditional$Freq.y

  matr <- merge(absolute, conditional, by = 1, all = T)
  matr <- matr[matr$Var1 != "0", ]
  matr[is.na(matr)] <- 0

  npidr <- matr$sum.y / matr$sum.x

  names(npidr) <- matr$Var1

  if(pool.type == "sum") {
    spool = apply(data, 1, sum)
  } else {
    spool = apply(data, 1, max)
  }

  npidr_tab <- npidr[as.character(spool)]

  npidr_tab[is.na(npidr_tab)] <- 0

  return(npidr_tab)
}