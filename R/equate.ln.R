equate.ln <- function(x, y, type = "linear", ident = 0, method = NA,
	w = -1, internal = TRUE, lts = FALSE, verbose = FALSE, ...) {

	xscale <- unique(x[, 1])
	yscale <- unique(y[, 1])
	xcount <- as.vector(tapply(x[, ncol(x)], x[, 1], sum))
	ycount <- as.vector(tapply(y[, ncol(x)], y[, 1], sum))
	type <- match.arg(tolower(type), c("mean", "linear"))
	method <- match.arg(tolower(method),
		c(NA, "nominal weights", "tucker", "levine", "chained",
		"braun/holland"))

	if(is.na(method)) {
		slope <- ifelse(type == "mean", 1,
			sd.freqtab(y)/sd.freqtab(x))
		intercept <- mean(y) - slope*mean(x)
	}
	else if(method == "chained") {
		if(type == "mean")
			slope1 <- slope2 <- 1
		else {
			slope1 <- sd.freqtab(y[, -2])/sd.freqtab(y[, -1])
			slope2 <- sd.freqtab(x[, -1])/sd.freqtab(x[, -2])
		}
		intercept <- mean.freqtab(y) -
			slope1*(slope2*mean.freqtab(x) + mean.freqtab(x[, -1]) -
			mean.freqtab(y[, -1]))
		slope <- slope1*slope2
	}
	else {
		stabs <- synthetic(x, y, w, method, internal, lts)
		stats <- stabs[[1]]
		if(!lts) {
			slope <- ifelse(type == "mean", 1, stats[4]/stats[3])
			intercept <- stats[2] - slope*stats[1]
		}
		else {
			slope <- ifelse(type == "mean", 1, stats[2]/stats[1])
			intercept <- mean(y) - slope*mean(x) +
				stats[2]*(mean.freqtab(x[, -1]) - mean.freqtab(y[, -1]))
		}
	}

	yx <- slope*xscale + intercept
 	yx <- (1 - ident)*yx + ident*xscale

	if(verbose) {
		out <- list(yx = yx)
		out$stats <- rbind(x = c(descript(x)), y = c(descript(y)),
			yx = c(descript(as.freqtab(cbind(yx, xcount)))))
		colnames(out$stats) <- c("mean", "sd", "skew", "kurt", "n")
		out$x <- x
		out$y <- y
		out$coefficients <- rbind(intercept, slope)[, 1]
		if(is.na(method))
			out$se <- se.ln(x, y)
		else {
			out$anchorstats <- rbind(descript(x[, -1]),
				descript(y[, -1]))
			rownames(out$anchorstats) <- c("xv", "yv")
			out$anchortab <- cbind(scale = unique(x[, 2]),
				xvcount = tapply(x[, 3], x[, 2], sum),
				yvcount = tapply(y[, 3], y[, 2], sum))
			if(method != "chained" & !lts)
				out <- c(out, stabs)
		}
	}
	else out <- yx

	return(out)
}
