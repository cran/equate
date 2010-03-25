max.freqtab <- function(x, ..., na.rm = FALSE) {
  max(x[as.logical(x[, ncol(x)]), 1])
}
