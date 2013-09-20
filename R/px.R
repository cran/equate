px <- function(x, y) {

	if(class(x) == "freqtab" & missing(y)) {
		colnames(x) <- NULL
		x[, 2] <- x[, 2]/sum(x[, 2])
		p <- .5*x[1, 2]
		for(i in 2:nrow(x))
			p[i] <- sum(x[1:(i - 1), 2]) + .5*x[i, 2]
	}
	else {
		x <- cbind(x)[, 1]
		xs <- floor(x + .5)
		yn <- sum(y[, 2])
		f <- sapply(xs, function(x)
			sum(y[y[, 1] <= x, 2])/yn)
		flow <- sapply(xs, function(x)
			sum(y[y[, 1] <= x - 1, 2])/yn)
		p <- flow + (x - (xs - .5))*(f - flow)
	}

	return(p)
}
