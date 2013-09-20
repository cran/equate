min.freqtab <- function(x, ..., na.rm = FALSE) {

	return(min(x[as.logical(x[, ncol(x)]), 1]))
}
