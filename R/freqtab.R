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

	if(missing(xscale))
		xscale <- min(xt):max(xt)
	if(is.null(vt)) {
		out <- cbind(xscale,
			as.vector(table(factor(xt, levels = xscale))))
	}
	else {
		if(missing(vscale))
			vscale <- min(vt):max(vt)
		temptab <- table(factor(vt, levels = vscale),
			factor(xt, levels = xscale))
		out <- cbind(rep(xscale, each = length(vscale)),
			rep(vscale, length(xscale)), as.vector(temptab))
	}
	out <- as.freqtab(out)

  return(out)
}
