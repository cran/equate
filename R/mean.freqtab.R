mean.freqtab <- function(x, ...) {
  sum(x[, 1] * x[, ncol(x)])/sum(x[, ncol(x)])
}
