#----------------------------------------------------------------
# Methods for presmoothing frequency distributions

#----------------------------------------------------------------
# Main presmoothing function

presmoothing <- function(x, smoothmethod = c("none",
	"average", "bump", "loglinear"), jmin,
	asfreqtab = TRUE, ...) {

	smoothmethod <- match.arg(smoothmethod)
	if(smoothmethod == "none")
		return(x)
	else if(smoothmethod == "average")
		return(freqavg(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if(smoothmethod == "bump")
		return(freqbump(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if(smoothmethod == "loglinear")
		return(loglinear(x, asfreqtab = asfreqtab, ...))
}

#----------------------------------------------------------------
# Loglinear smoothing
# Support coming for comparing a vector of degrees

loglinear <- function(x, scorefun, degree = 4, xdegree = 1,
	raw = TRUE, asfreqtab = TRUE, verbose = FALSE,
	compare = FALSE, stepup = compare, showWarnings = FALSE,
	...) {

	nc <- ncol(x)
	if(missing(scorefun)) {
		scorefun <- NULL
		for(i in 1:(nc - 1)) {
			tempfun <- poly(x[, i],
				degree = degree, raw = raw)
			colnames(tempfun) <- paste(c("x", "v")[i],
				1:degree, sep = "")
			scorefun <- cbind(scorefun, tempfun)
		}
		if(nc == 3 & xdegree > 0) {
			tempfun <- poly(x[, 1]*x[, 2],
				degree = xdegree, raw = raw)
			colnames(tempfun) <- paste("xv",
				1:xdegree, sep = "")
			scorefun <- cbind(scorefun, tempfun)
		}
	}
	else if(nrow(scorefun) != nrow(x))
		stop("'scorefun' must contain the same ",
			"number of rows as 'x'")
	scorefun <- data.frame(f = x[, nc], scorefun)

	if(ncol(scorefun) < 3 & (stepup | compare))
		stop("to run multiple models, 'scorefun' must",
			" include multiple variables")

	if(stepup | compare) {
		if(showWarnings)
			out <- lapply(2:ncol(scorefun), function(i)
				glm(scorefun[, 1:i], family = poisson))
		else
			suppressWarnings(out <- lapply(2:ncol(scorefun),
				function(i) glm(scorefun[, 1:i],
					family = poisson)))
		names(out) <- colnames(scorefun)[-1]
	}
	else if(showWarnings)
		out <- glm(scorefun, family = poisson)
	else
		suppressWarnings(out <- glm(scorefun, family = poisson))

	if(compare) {
		nm <- length(out)
		resdf <- as.numeric(lapply(out, function(x) x$df.residual))
		resdev <- as.numeric(lapply(out, function(x) x$deviance))
		aic <- as.numeric(lapply(out, AIC))
		bic <- as.numeric(lapply(out, BIC))
		tab <- data.frame(resdf, resdev, aic, bic,
			c(NA, -diff(resdf)), c(NA, 
			-diff(resdev)))
		vars <- lapply(out, function(x) paste(deparse(formula(x)), 
			collapse = "\n"))
		dimnames(tab) <- list(1:nm, c("Resid. Df", "Resid. Dev", 
			"AIC", "BIC", "Df", "Deviance"))
		tab <- stat.anova(table = tab, test = "Chisq", scale = 1, 
			df.scale = Inf,
			n = length(out[[order(resdf)[1]]]$residuals))
		return(structure(tab,
			heading = c("Analysis of Deviance Table\n",
			paste("Model ", format(1:nm), ": ", vars, 
			sep = "", collapse = "\n")), class = c("anova", 
			"data.frame")))
	}
	else if(verbose)
		return(out)
	else if(stepup)
		return(data.frame(lapply(out, fitted)))
	else if(asfreqtab)
		return(as.freqtab(cbind(x[, -ncol(x)],
			out$fitted)))
	else
		return(out$fitted)
}

#----------------------------------------------------------------
# Frequency adjustment
# Bump frequencies upward by a small amount

freqbump <- function(x, jmin = 1e-6, asfreqtab = FALSE, ...) {

	nc <- ncol(x)
	f <- x[, nc]/sum(x[, nc])
	out <- (f + jmin)/(1 + (max(x[, 1]) + 1)*jmin)
	out <- ((f + jmin)/sum(f + jmin))*sum(x[, nc])
	
	if(asfreqtab)
		return(as.freqtab(cbind(x[, -nc], out)))
	else
		return(out)
}

#----------------------------------------------------------------
# Frequency averaging

freqavg <- function(x, jmin = 1, asfreqtab = FALSE, ...) {
	
	xtab <- x
	nc <- ncol(xtab)
	if(nc > 2)
		stop("frequency averaging only supported for univariate 'x'")
		
	x <- cbind(x, 0, 0, 0, 0, 0)
	ks <- 1
	while(ks <= nrow(x) & x[ks, 2] < jmin)
		ks <- ks + 1
	x[1:ks, 3] <- 1
	x[1:ks, 4] <- ks

	lls <- nrow(x)
	while(lls >= 0 & x[lls, 2] < jmin)
		lls <- lls - 1
	x[lls:nrow(x), 3] <- lls
	x[lls:nrow(x), 4] <- nrow(x)

	ss <- ks + 1
	tts <- lls - 1
	for(j in ss:tts) {
		if(x[j, 2] < jmin) {
			lls <- j
			ks <- j
			while(lls >= 1 & x[lls, 2] < jmin)
				lls <- lls - 1
			while(ks <= nrow(x) & x[ks, 2] < jmin)
				ks <- ks + 1
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(p in 1:(nrow(x) - 1)) {
		if(x[p, 4] > 0 & x[p, 4] == x[p + 1, 3]) {
			if(x[p, 3] > 0)
				lls <- x[p, 3]
			if(x[p + 1, 4] > 0)
				ks <- x[p + 1, 4]
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(j in 1:nrow(x)) {
		lls <- x[j, 3]
		if(lls == 0)
			lls <- j
		ks <- x[j, 4]
		if(ks == 0)
			ks <- j
		sumit <- 0
		sumit <- sumit + sum(x[lls:ks, 2])
		for(i in lls:ks) {
			x[i, 5] <- sumit
			x[i, 6] <- x[j, 4] - x[j, 3] + 1
			x[i, 7] <- x[i, 5]/x[i, 6]
		}
		j <- j + x[j, 4] - x[j, 3]
	}

	colnames(x)[c(1, 2, 6, 7)] <-
		c("score", "count", "b", "acount")

	if(asfreqtab)
		return(as.freqtab(cbind(xtab[, -nc],
			x[, 7])))
	else
		return(x[, 7])
}

#----------------------------------------------------------------
