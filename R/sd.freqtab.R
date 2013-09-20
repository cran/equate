sd.freqtab <- function(x) {

	return(sqrt(cov.freqtab(x[, c(1, ncol(x))])))
}
