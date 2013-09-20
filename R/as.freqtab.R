as.freqtab <- function(x) {

	if(ncol(x) == 3)
		colnames(x) <- c("x", "v", "count")
	else if(ncol(x) == 2)
		colnames(x) <- c("x", "count")
	else
		stop("incorrect number of columns")
	rownames(x) <- NULL
	class(x) <- "freqtab"

	return(x)
}
