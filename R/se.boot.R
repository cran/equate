se.boot <- function(x, y, xn = sum(x[, ncol(x)]),
	yn = sum(y[, ncol(y)]), reps = 100, eqfun,
	returnboots = FALSE, ...) {

	nc <- ncol(x)
	xscale <- unique(x[, 1])
	yscale <- unique(y[, 1])
	xprob <- x[, nc]/sum(x[, nc])
	yprob <- y[, nc]/sum(y[, nc])
	eqfun <- match.fun(eqfun)
	eqmat <- matrix(nrow = length(xscale), ncol = reps)
	if(nc == 3) {
		index <- 1:nrow(x)
		vscale <- unique(x[, 2])
		for(i in 1:reps) {
			xi <- sample(index, xn, replace = TRUE, prob = xprob)
			xtemp <- freqtab(x[xi, 1], x[xi, 2], xscale = xscale,
				vscale = vscale)
			yi <- sample(index, yn, replace = TRUE, prob = yprob)
			ytemp <- freqtab(y[yi, 1], y[yi, 2], xscale = yscale,
				vscale = vscale)
			eqmat[, i] <- eqfun(xtemp, ytemp, ...)
		}
	}
	else {
		for(i in 1:reps) {
			xtemp <- freqtab(sample(x[, 1], xn,
				replace = TRUE, prob = xprob), xscale = xscale)
			ytemp <- freqtab(sample(y[, 1], yn,
				replace = TRUE, prob = yprob), xscale = yscale)
			eqmat[, i] <- eqfun(xtemp, ytemp, ...)
		}
	}

	if(returnboots)
		return(eqmat)
	else
		return(apply(eqmat, 1, sd))
}
