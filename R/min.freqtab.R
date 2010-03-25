min.freqtab <- function(x, ..., na.rm = FALSE) {
  min(x[as.logical(x[, ncol(x)]), 1])
}
