plot.equate <- function(..., add = FALSE, addident = TRUE,
	xpoints, ypoints, xlab = "Raw Scale", ylab = "Equated Scale",
	main = "Equating Plot", col = rainbow(length(x)),
	lty = 1, pchcol = "gray", addlegend = TRUE,
	legendtext, morepar = NULL) {
	
	x <- list(...)
	if(class(x) == "equate")
		x <- list(x)
	nx <- length(x)
	xscale <- unique(x[[1]]$x[, 1])
	yscale <- unique(x[[1]]$y[, 1])
	
	if(!add) {
		plot(c(min(xscale), max(xscale)),
			c(min(yscale), max(yscale)), xlab = xlab,
			ylab = ylab, main = main, type = "n", unlist(morepar))
		if(!missing(xpoints) & !missing(ypoints))
			points(xpoints, ypoints, col = pchcol)
	}
	
	if(addident)
		lines(xscale, xscale)
	
	lty <- rep(lty, length = nx)
	col <- rep(col, length = nx)	
	for(i in 1:nx) {
		lines(x[[i]]$concordance, col = col[i],
			lty = lty[i], unlist(morepar))
	}
			
	if(addlegend) {
		if(missing(legendtext)) {
			legendtext <- c("mean", "linear", "circle", "equip")
			legendtext <- lapply(x, function(y)
				legendtext[charmatch(substr(y$type, 1, 1),
				legendtext)])
			if(x[[1]]$design == "nonequivalent groups") {
				methods <- c("nW", "chain", "b/H", "tucker",
					"levine", "fE")
				methods <- lapply(x, function(y)
					methods[charmatch(substr(y$method, 1, 1),
					methods)])
				legendtext <- paste(legendtext,
					methods, sep = ": ")				
			}
			legendtext <- gsub("\\b(\\w)", "\\U\\1",
				legendtext, perl = TRUE)
		}
		if(addident) {
			legendtext <- c("Identity", legendtext)
			lty = c(1, lty)
			col = c(1, col)
		}
		legend("bottomright", legend = legendtext,
			lty = lty, col = col, bty = "n")
	}
}