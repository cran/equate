mean.freqtab <- function(x, ...) {

	return(sum(x[, 1]*x[, ncol(x)])/sum(x[, ncol(x)]))
}
