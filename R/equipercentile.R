#----------------------------------------------------------------
# Equipercentile equating functions

#----------------------------------------------------------------
# The main function

equipercentile <- function(x, y, method = "none",
	ws = -1, smoothmethod = c("none", "bump",
	"average", "loglinear"), lowp = c(min(x[, 1]),
	min(y[, 1])), highp = c(max(x[, 1]), max(y[, 1])),
	verbose = FALSE, ...) {

	smoothmethod <- match.arg(smoothmethod)
	if(missing(y))
		y <- mfreqtab(x, 2)
	if(ncol(y) < ncol(x)) {
		xtab <- presmoothing(x, smoothmethod, ...)
		ytab <- mfreqtab(xtab, 2)
		xtab <- mfreqtab(xtab)
	}
	else {
		xtab <- presmoothing(x, smoothmethod, ...)
		ytab <- presmoothing(y, smoothmethod, ...)
	}

	if(method == "none")
		yxs <- equip(xtab, ytab, highp[2])
	else if(method == "frequency estimation") {
		stabs <- synthetic(xtab, ytab, ws, method)
		yxs <- equip(mfreqtab(stabs$xsynthetic),
			mfreqtab(stabs$ysynthetic), highp[2])
	}
	else if(method == "chained") {
		vx <- equip(mfreqtab(xtab),
			mfreqtab(xtab, 2))$yx
		pvyx <- px(vx, mfreqtab(ytab, 2))
		yxs <- equip(pvyx, mfreqtab(ytab))
	}
	yx <- yxs$yx
	
	if(!verbose)
		out <- yx
	else {
		out <- c(list(x = x, y = y, concordance = data.frame(scale =
			unique(x[, 1]), yx = yx), points = data.frame(low = lowp,
			high = highp, row.names = c("x", "y")), smoothmethod =
			smoothmethod), list(...)[names(list(...)) %in% c("jmin",
			"degree", "xdegree", "scorefun")])
		if(method == "frequency estimation")
			out <- c(out, stabs)
		if(smoothmethod != "none") {
			out$xsmooth <- xtab
			out$ysmooth <- ytab
		}
	}
	return(out)
}

#----------------------------------------------------------------
# Equipercentile conversion of x to y
# Does this work for scales with negatives? Yes
# Does it work for shrunken or stretched scales? Yes, now.
# Previously, it is assumed that scores were in 1 point increments

equip <- function(x, y, ky = max(y[, 1])) {

	yscale <- unique(y[, 1])
	yn <- sum(y[, 2])
	if(!is.freqtab(x)) {
		prank <- sort(unique(x))
		xscale <- yscale
		xn <- 0
	}
	else {
		prank <- px(x)
		xscale <- unique(x[, 1])
		xn <- sum(x[, 2])
		se <- vector(length = length(prank))
	}
	yinc <- round(diff(yscale)[1], 8)
	if(any(round(diff(yscale), 8) != yinc))
		stop("'y' scale must be equal-interval")
	hinc <- yinc/2
	yx <- vector()
	fy <- fx(y)
	sn <- nrow(y)
	Ly <- min(yscale)
	xnone <- prank == 0
	xone <- prank == 1
	xbot <- sum(xnone) + 1
	xtop <- sum(!xone)
	xyone <- which(xscale[xone] > (ky + hinc)) + xtop

	yx[xnone] <- Ly - hinc
	yx[xone] <- ky + hinc
	yx[xyone] <- xscale[xyone]
	for(j in xbot:xtop) {
		yu <- 1
		while(fy[yu] <= prank[j])
			yu <- yu + 1
		yu2 <- ifelse(yu == 1, 0, fy[yu - 1])
		g0 <- fy[yu] - yu2
		yx[j] <- y[yu, 1] - hinc + ((prank[j] - yu2)/g0)*yinc
		if(g0 > 0 & xn)
			se[j] <- eqse(prank[j], g0, yu2, xn, yn)
	}
	if(any(y[, 2] == 0)) {
		xbot <- xbot + sum(prank[!xnone] <= min(fy))
		for(i in xbot:xtop) {
			yl <- sn
			while(fy[yl] >= prank[i])
				yl <- yl - 1
			yl2 <- fy[yl + 1]
			yxtemp <- y[yl, 1] + hinc + ((prank[i] - fy[yl])/
				(yl2 - fy[yl]))*yinc
			yx[i] <- mean(c(yx[i], yxtemp))
		}
	}

	if(xn) {
		out <- list(yx = yx)
		out$se <- se
	}
	else
		out <- list(yx = yx[match(x, prank)])

	return(out)
}

#----------------------------------------------------------------
