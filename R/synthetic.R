#----------------------------------------------------------------
# Functions for obtaining frequencies and parameters
# for the synthetic distribution

#----------------------------------------------------------------
# Main function

synthetic <- function(x, y, ws = .5, method, internal = TRUE,
	lts = FALSE) {

	xscale <- unique(x[, 1])
	yscale <- unique(y[, 1])
	vscale <- unique(x[, 2])
	if(ws == -1)
		ws <- sum(x[, 3])/sum(x[, 3], y[, 3])
	yws <- 1 - ws

	if(method == "frequency estimation" |
		method == "braun/holland") {
		synth <- syntheticequip(x, y, ws)
		fstab <- as.freqtab(cbind(x[, -3],
			synth$fs))
		gstab <- as.freqtab(cbind(y[, -3],
			synth$gs))
		if(method == "braun/holland") {
			msx <- mean(mfreqtab(fstab))
			msy <- mean(mfreqtab(gstab))
			sdsx <- sd.freqtab(mfreqtab(fstab))
			sdsy <- sd.freqtab(mfreqtab(gstab))
		}
	}
	else {
		mx <- mean(x)
		varx <- cov.freqtab(x[, -2])
		mxv <- mean.freqtab(x[, -1])
		varxv <- cov.freqtab(x[, -1])
		my <- mean(y)
		vary <- cov.freqtab(y[, -2])
		myv <- mean.freqtab(y[, -1])
		varyv <- cov.freqtab(y[, -1])

		covxv <- cov.freqtab(x)
		covyv <- cov.freqtab(y)

		if(method == "nominal weights")
			g1 <- g2 <- max(xscale)/max(vscale)
		else if(method == "tucker") {
			g1 <- covxv/varxv
			g2 <- covyv/varyv
		}
		else if(method == "levine" & internal) {
			g1 <- varx/covxv
			g2 <- vary/covyv
		}
		else if(method == "levine" & !internal) {
			g1 <- (varx + covxv)/(varxv + covxv)
			g2 <- (vary + covyv)/(varyv + covyv)
		}
		if(!lts) {
			msx <- mx - (yws*g1*(mxv - myv))
			msy <- my + (ws*g2*(mxv - myv))
			sdsx <- sqrt(varx -
				(yws*(g1^2)*(varxv - varyv)) +
				(ws*yws*(g1^2)*(mxv - myv)^2))
			sdsy <- sqrt(vary +
				(ws*(g2^2)*(varxv - varyv)) +
				(ws*yws*(g2^2)*(mxv - myv)^2))
		}
	}
	if(lts)
		out <- list(gamma = c(g1, g2))
	else {
		out <- list(ws = ws)
		if(method == "frequency estimation" |
			method == "braun/holland")
			out <- c(out, list(xsynthetic = fstab,
				ysynthetic = gstab))
		if(method != "frequency estimation") {
			out$synthstats <- data.frame(mean = c(msx, msy),
				sd = c(sdsx, sdsy))
			rownames(out$synthstats) <-
				c("xsynthetic", "ysynthetic")
		}
	}
	return(out)
}

#----------------------------------------------------------------
# Frequency estimation

syntheticequip <- function(x, y, ws) {

	h1 <- mfreqtab(x, 2)[, 2]
	f1xv <- x[, 3]/h1
	f1xv[is.na(f1xv)] <- 0
	h2 <- mfreqtab(y, 2)[, 2]
	g2yv <- y[, 3]/h2
	g2yv[is.na(g2yv)] <- 0
	f2 <- f1xv*h2
	g1 <- g2yv*h1
	fs <- ws*x[, 3] + (1 - ws)*f2
	gs <- ws*g1 + (1 - ws)*y[, 3]

	return(list(fs = fs, gs = gs))
}

#----------------------------------------------------------------
