sd.freqtab <- function(x) {
  sqrt(cov.freqtab(x[, c(1, ncol(x))]))
}
