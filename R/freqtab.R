#----------------------------------------------------------------
# Functions for creating, manipulating, summarizing,
# and plotting frequency tables

#----------------------------------------------------------------
# Main function for creating a frequency table

freqtab <- function(x, v, xitems = 1:ncol(x),
	vitems = 1:ncol(v), xscale, vscale, na.rm = TRUE) {

	x <- cbind(x)
	
	if(ncol(x) > 1) {
		if(is.factor(xitems))
			xitems <- as.character(xitems)
		xt <- apply(x[, xitems], 1, sum, na.rm = na.rm)		
	}
	else
		xt <- x

	if(!missing(v)) {
		if(ncol(cbind(v)) > 1) {
			if(is.factor(vitems))
				vitems <- as.character(vitems)
			vt <- apply(v[, vitems], 1, sum, na.rm = na.rm)
		}
		else
			vt <- v
	}
	else if(!missing(vitems)) {
		if(is.factor(vitems))
			vitems <- as.character(vitems)
		vt <- apply(x[, vitems], 1, sum, na.rm = na.rm)
	}
	else
		vt <- NULL

	xtc <- xt[complete.cases(xt, vt)]
	vtc <- vt[complete.cases(xt, vt)]
	
	if(missing(xscale))
		xscale <- min(xtc):max(xtc)
	if(is.null(vt)) {
		out <- cbind(xscale,
			as.vector(table(factor(xtc, levels = xscale))))
	}
	else {
		if(missing(vscale))
			vscale <- min(vtc):max(vtc)
		temptab <- table(factor(vtc, levels = vscale),
			factor(xtc, levels = xscale))
		out <- cbind(rep(xscale, each = length(vscale)),
			rep(vscale, length(xscale)), as.vector(temptab))
	}
	out <- as.freqtab(out)

  return(out)
}

#----------------------------------------------------------------
# Assign freqtab class

as.freqtab <- function(x) {

	x <- as.data.frame(unclass(x))
	if(ncol(x) == 3)
		colnames(x) <- c("total", "anchor", "observed")
	else if(ncol(x) == 2)
		colnames(x) <- c("total", "observed")
	else
		stop("incorrect number of columns")
	rownames(x) <- NULL
	class(x) <- c("freqtab", "data.frame")

	return(x)
}

#----------------------------------------------------------------
# Test for freqtab class

is.freqtab <- function(x) {
	
	return(class(x)[1] == "freqtab")
}

#----------------------------------------------------------------
# Convert bivariate frequency table to univariate by
# summing over one margin

mfreqtab <- function(x, margin = 1) {

	if(!is.freqtab(x))
		stop("class of 'x' must be 'freqtab'")
	return(as.freqtab(cbind(unique(x[, margin]),
		tapply(x[, ncol(x)], x[, margin], sum))))
}

#----------------------------------------------------------------
# Plot method

plot.freqtab <- function(x, y = NULL, xcol = 1,
	ycol, pch = 16, ylty = 1, xlab = "Total Test",
	addlegend = !missing(y), legendtext, ...) {
	
	if(ncol(x) == 2)
		ufreqplot(x, y, xcol, ycol, ylty, xlab,
			horiz = FALSE, addlegend = addlegend,
			legendtext = legendtext, ...)
	else if(ncol(x) == 3)
		bfreqplot(x, y, xcol, ycol, pch, ylty, xlab,
			addlegend = addlegend,
			legendtext = legendtext, ...)
	else stop("'x' must be either univariate or bivariate")
}

#----------------------------------------------------------------
# Internal univariate plot

ufreqplot <- function(x, y = NULL, xcol = 1, ycol,
	ylty = 1, xlab = "Total Test", ylab = "Count",
	horiz = FALSE, addlegend = FALSE,
	legendtext, ...) {

	if(!is.null(y)) {
		if(is.freqtab(y))
			y <- cbind(y[, 2])
		else
			y <- cbind(y)
		if(missing(ycol))
			ycol <- rainbow(ncol(y))
	}
	
	if(horiz) {
		plot.default(round(range(0, x[, 2], y)),
			range(x[, 1]), type = "n", xlab = xlab,
			ylab = ylab, ...)
		segments(rep(0, nrow(x)), x[, 1], x[, 2],
			col = xcol)
		if(!is.null(y))
			matlines(y, x[, 1], col = ycol,
				lty = ylty, ...)
	}
	else {
		plot.default(range(x[, 1]),
			round(range(0, x[, 2], y)), type = "n",
			xlab = xlab, ylab = ylab, ...)
		segments(x[, 1], y0 = rep(0, nrow(x)),
			y1 = x[, 2], col = xcol)
		if(!is.null(y))
			matlines(x[, 1], y, col = ycol,
				lty = ylty, ...)
	}
	
	if(addlegend & !is.null(y)) {
		if(missing(legendtext))
			legendtext <- if(is.null(colnames(y)))
				1:ncol(y) else colnames(y)
		legend("topright", legend = legendtext,
			lty = ylty, col = ycol, bty = "n")
	}
}

#----------------------------------------------------------------
# Internal bivariate plot

bfreqplot <- function(x, y = NULL, xcol = 1,
	ycol, pch = 16, ylty = 1, xlab = "Total Test",
	ylab = "Anchor Test", addlegend = FALSE,
	legendtext, ...) {

	if(!is.null(y)) {
		if(is.freqtab(y))
			y <- cbind(y[, 3])
		if(missing(ycol))
			ycol <- rainbow(ncol(y))
		ytab <- apply(y, 2, function(z)
			tapply(z, x[, 1], sum))
		yvtab <- apply(y, 2, function(z)
			tapply(z, x[, 2], sum))
	}
	else ytab <- yvtab <- NULL
	
	reset.par <- par(no.readonly = TRUE)
	nf <- layout(matrix(c(2, 4, 1, 3), 2, 2,
		byrow = TRUE), c(3, 1), c(1, 3), TRUE)
	par(mar = c(4, 4, 1, 1))
	plot(range(x[, 1]), range(x[, 2]), type = "n",
		xlab = xlab, ylab = ylab, ...)
	points.freqtab(x, xcol = xcol, pch = pch)

	par(mar = c(0, 4, 1, 1))
	xtab <- mfreqtab(x)
	ufreqplot(xtab, ytab, xcol, ycol, ylty,
		xlab = "", ylab = "", xaxt = "n", bty = "n")

	par(mar = c(4, 0, 1, 1))
	xvtab <- mfreqtab(x, 2)
	ufreqplot(xvtab, yvtab, xcol, ycol, ylty,
		xlab = "", ylab = "", yaxt = "n", bty = "n",
		horiz = TRUE)

	if(addlegend & !is.null(y)) {
		par(mar = c(0, 0, 0, 0))
		plot(0, 0, type = "n", bty = "n", xaxt = "n",
			yaxt = "n")
		if(missing(legendtext))
			legendtext <- if(is.null(colnames(y)))
				1:ncol(y) else colnames(y)
		legend("bottomleft", legend = legendtext,
			lty = ylty, col = ycol, bty = "n")
	}
	
	par(reset.par)
}

#----------------------------------------------------------------
# Points method

points.freqtab <- function(x, xcol = 1, pch = 16, ...) {
	
	if(ncol(x) < 3)
		stop("'x' must be a bivariate frequency table")
		
	index <- as.logical(x[, 3])
	xpoints <- x[index, 1]
	vpoints <- x[index, 2]
	if(sd(x[index, 3]) > 0)
		dens <- pmax(0, pmin(255,
			scale(x[index, 3])*50 + 100))
	else dens <- rep(150, sum(index))
	rgbcol <- col2rgb(xcol)
	ptcol <- rgb(rgbcol[1], rgbcol[2], rgbcol[3],
		dens, maxColorValue = 255)
	points(xpoints, vpoints, col = ptcol,
		pch = pch, ...)	
}

#----------------------------------------------------------------
# Summary method

summary.freqtab <- function(object, ...) {
	
	nc <- ncol(object) - 1
	out <- NULL
	for(i in 1:nc)
		out <- rbind(out,
			descript(object[, c(i, nc + 1)]))
	rownames(out)[1:nc] <- c("total", "anchor")[1:nc]
	return(out)
}

#----------------------------------------------------------------
# Internal descriptives function
# For univariate frequency table

descript <- function(x) {
	if(is.freqtab(x))
		return(data.frame(
			mean = mean.freqtab(x),
			sd = sd.freqtab(x),
			skew = skew.freqtab(x),
			kurt = kurt.freqtab(x),
			min = min.freqtab(x),
			max = max.freqtab(x),
			n = sum(x[, 2])))
}

#----------------------------------------------------------------
# Mean

mean.freqtab <- function(x, ...) {

	return(sum(x[, 1]*x[, ncol(x)])/sum(x[, ncol(x)]))
}

#----------------------------------------------------------------
# Standard deviation

sd.freqtab <- function(x) {

	return(sqrt(cov.freqtab(x[, c(1, ncol(x))])))
}

#----------------------------------------------------------------
# Covariance

cov.freqtab <- function(x) {

	nc <- ncol(x)
	i <- x[, nc] > 0
	xc <- x[i, 1] - mean.freqtab(x)
	vc <- x[i, nc - 1] - mean.freqtab(x[, (nc - 1):nc])

	return(sum(xc*vc*x[i, nc])/(sum(x[, nc]) - 1))
}

#----------------------------------------------------------------
# Correlation

cor.freqtab <- function(x) {

	return(cov.freqtab(x)/(sd.freqtab(x[, -2])*sd.freqtab(x[, -1])))
}

#----------------------------------------------------------------
# Minimum

min.freqtab <- function(x, ..., na.rm = FALSE) {

	return(min(x[as.logical(x[, ncol(x)]), 1]))
}

#----------------------------------------------------------------
# Maximum

max.freqtab <- function(x, ..., na.rm = FALSE) {

	return(max(x[as.logical(x[, ncol(x)]), 1]))
}

#----------------------------------------------------------------
# Skewness

skew.freqtab <- function(x) {

	nc <- ncol(x)
	i <- x[, nc] > 0

	return(sum(((x[i, 1] - mean.freqtab(x))^3)*x[i, nc])/
		(sum(x[i, nc]) - 1)/cov.freqtab(x[, c(1, nc)])^1.5)
}

#----------------------------------------------------------------
# Kurtosis

kurt.freqtab <- function(x) {

	nc <- ncol(x)
	i <- x[, nc] > 0

	return(sum(((x[i, 1] - mean.freqtab(x))^4)*x[i, nc])/
		(sum(x[i, nc]) - 1)/cov.freqtab(x[, c(1, nc)])^2)
}

#----------------------------------------------------------------
# Cumulative frequency

fx <- function(x) {

	if(!is.null(dim(x)))
		x <- x[, ncol(x)]
	x <- x/sum(x)

	f <- numeric(length(x))
	for(i in 1:length(x))
		f[i] <- sum(x[1:i])

	return(f)
}

#----------------------------------------------------------------
# Percentile ranks

px <- function(x, y) {

	if(is.freqtab(x) & missing(y)) {
		x[, 2] <- x[, 2]/sum(x[, 2])
		p <- .5*x[1, 2]
		for(i in 2:nrow(x))
			p[i] <- sum(x[1:(i - 1), 2]) + .5*x[i, 2]
	}
	else if(is.freqtab(y)) {
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

#----------------------------------------------------------------
