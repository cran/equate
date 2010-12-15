kurt.freqtab <- function(x) {
  nc <- ncol(x)
  i <- x[, nc] > 0
  sum(((x[i, 1] - mean.freqtab(x))^4) * x[i, nc])/
    (sum(x[i, nc]) - 1)/cov.freqtab(x[, c(1, nc)])^2
}
