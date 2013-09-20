se.ln <- function(x, y) {

	nx <- sum(x[, 2])
	ny <- sum(y[, 2])
	mx <- mean(x)
	sdx <- sd.freqtab(x)
	vary <- cov.freqtab(y)
	skewterm <- skew.freqtab(x)/nx + skew.freqtab(y)/ny
	kurtterm <- (kurt.freqtab(x) - 1)/(4*nx) +
		(kurt.freqtab(y) - 1)/(4*ny)
	se <- vector()
	for(i in 1:nrow(x)) {
		xterm <- (x[i, 1] - mx)/sdx
		se[i] <- vary*(1/nx + 1/ny + skewterm*xterm +
			kurtterm*xterm^2)
	}

	return(se)
}
