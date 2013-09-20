equate <- function(x, y, type, method = NA, name = NULL,
	ident = 0, bootse = FALSE, ...) {

	if(class(y) == "equate") {
		if(y$type == "equipercentile") {
			if(is.na(y$method) | y$method == "chained") {
				xtab <- as.freqtab(cbind(y$x[, 1],
					y$x[, ncol(y$x)]))
				ytab <- as.freqtab(cbind(y$y[, 1],
					y$y[, ncol(y$y)]))
			}
			else {
				xtab <- y$xsynthtab
				ytab <- y$ysynthtab
			}
			if(is.na(y$method) | y$method == "frequency estimation")
				p <- px(x, xtab)
			else {
				xx <- px(x, xtab)
				xv <- equipercentile(xx,
					y$anchortab[, c(1, ncol(y$anchortab) - 1)])$yx
				p <- px(xv,
					y$anchortab[, c(1, ncol(y$anchortab))])
			}
			out <- equipercentile(p, ytab)$yx
		}
		else if(y$type == "circle-arc") {
			index <- x >= y$points[1] & x <= y$points[5]
			out <- lin <- y$coef[2]*x + y$coef[1]
			if(y$points[4] < y$points[3])
				out[index] <- lin[index] + (y$coef[4] -
					sqrt((y$coef[5]^2) - (lin[index] - y$coef[3])^2))
			else if(y$points[4] > y$points[3])
				out[index] <- lin[index] + (y$coef[4] +
					sqrt((y$coef[5]^2) - (lin[index] - y$coef[3])^2))
		}
		else out <- y$coef[2]*x + y$coef[1]
		names(out) <- NULL

		return((1 - y$ident)*out + y$ident*x)
	}

	type <- match.arg(tolower(type),
		c("mean", "linear", "circle-arc", "equipercentile"))
	eqfun <- match.fun(switch(type, equipercentile = "equate.eq",
		"circle-arc" = "equate.ca", "equate.ln"))
	eqout <- eqfun(x, y, type = type, method = method,
		ident = ident, verbose = TRUE, ...)

	out <- list(name = name, type = type,
		method = match.arg(tolower(method),
			c(NA, "nominal weights", "tucker", "levine",
			"frequency estimation", "chained", "braun/holland")))
	out$design <- ifelse(is.na(method),
		"equivalent groups", "nonequivalent groups")
	out$ident <- ident
	out$concordance <- cbind(scale = unique(eqout$x[, 1]),
		yx = eqout$yx)
	if(bootse) {
		out$bootsee <- se.boot(x = x, y = y, eqfun = eqfun,
			type = type, method = method, ident = ident, ...)
	}
	out <- c(out, eqout)

	class(out) <- "equate"
	return(out)
}
