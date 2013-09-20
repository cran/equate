plot.freqtab <- function(x, y, col1 = "gray27",
	col2 = "gray74", type = "h", pch = 16, yoffset = .3, ...) {

	if(length(type) == 1)
		type[2] <- type
	
	if(ncol(x) == 2) {
		if(missing(y))
			plot.default(x, type = type[1], col = col1, ...)
		else {
			if(ncol(cbind(y)) == 1)
				y <- cbind(x[, 1], y)
			ylim <- round(c(min(min(x[, 2]), min(y[, 2])),
				max(max(x[, 2]), max(y[, 2]))))
			plot.default(x, type = type[1], ylim = ylim,
				col = col1, ...)
			points(y[, 1] + ifelse(type[2] == "h", yoffset, 0),
				y[, 2], type = type[2],
				col = col2, ...)
		}
	}
	else {
		reset.par <- par(no.readonly = TRUE)
		X <- c(min(x[, 1]), max(x[, 1]))
		V <- c(min(x[, 2]), max(x[, 2]))
		xcount <- as.vector(tapply(x[, 3], x[, 1], sum))
		vcount <- as.vector(tapply(x[, 3], x[, 2], sum))
		index <- as.logical(x[, 3])
		xpoints <- x[index, 1]
		vpoints <- x[index, 2]
		dens <- pmax(0, pmin(255, scale(x[index, 3]) * 50 + 100))
		ptcol <- rgb(50, 50, 50, dens, maxColorValue = 255)
		nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),
			c(3, 1), c(1, 3), TRUE)
		par(mar = c(4, 4, 1, 1))
		plot(X, V, type = "n", ...)
		points(xpoints, vpoints, col = ptcol, pch = pch, ...)
		par(mar = c(0, 3, 1, 1))
		barplot(xcount, axes = FALSE, space = 0, col = col1)
		par(mar = c(3, 0, 1, 1))
		barplot(vcount, axes = FALSE, space = 0, horiz = TRUE,
			col = col1)
		par(reset.par)
	}
}
