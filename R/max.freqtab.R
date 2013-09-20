max.freqtab <- function(x, ..., na.rm = FALSE) {

	return(max(x[as.logical(x[, ncol(x)]), 1]))
}
