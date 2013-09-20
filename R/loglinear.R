loglinear <- function(x, scorefun, degree, raw = TRUE,
	verbose = FALSE, compare = FALSE, stepup = compare,
	showWarnings = FALSE, ...) {

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
	}
	else if(nrow(scorefun) != nrow(x))
		stop("'scorefun' must contain the same number of rows as 'x'")
	scorefun <- data.frame(f = x[, nc], scorefun)

	if(ncol(scorefun) < 3 & (stepup | compare))
		stop("to run multiple models, 'scorefun' must include multiple variables")

	if(stepup | compare) {
		if(showWarnings)
			out <- lapply(2:ncol(scorefun), function(i)
				glm(scorefun[, 1:i], family = poisson, ...))
		else
			suppressWarnings(out <- lapply(2:ncol(scorefun), function(i)
				glm(scorefun[, 1:i], family = poisson, ...)))
		names(out) <- colnames(scorefun)[-1]
	}
	else if(showWarnings)
		out <- glm(scorefun, family = poisson, ...)
	else
		suppressWarnings(out <- glm(scorefun, family = poisson, ...))

	if(compare) {
		nm <- length(out)
		resdf <- as.numeric(lapply(out, function(x) x$df.residual))
		resdev <- as.numeric(lapply(out, function(x) x$deviance))
		tab <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
			-diff(resdev)))
		vars <- lapply(out, function(x) paste(deparse(formula(x)), 
			collapse = "\n"))
		dimnames(tab) <- list(1:nm, c("Resid. Df", "Resid. Dev", 
			"Df", "Deviance"))
		bigm <- out[[order(resdf)[1]]]
		tab <- stat.anova(table = tab, test = "Chisq", scale = 1, 
			df.scale = Inf, n = length(bigm$residuals))
		structure(tab, heading = c("Analysis of Deviance Table\n",
			paste("Model ", format(1:nm), ": ", vars, 
			sep = "", collapse = "\n")), class = c("anova", 
			"data.frame"))
	}
	else if(verbose)
		return(out)
	else if(stepup) {
		return(data.frame(lapply(out, fitted)))
	}
	else
		return(out$fitted)
}