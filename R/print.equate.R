print.equate <- function(x, ...) {

	up <- function(x)
		gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)

	design <- up(x$design)
	if(is.na(x$method)) {
		method <- up(x$type)
		stats <- x$stats
	}
	else {
		method <- paste(up(x$method), up(x$type))
		stats <- rbind(x$stats, x$anchorstats)
	}

	if(is.null(x$name))
		cat("\n", method, " Equating: ", design, sep = "")
	else
		cat(x$name)

	cat("\n\nSummary Statistics:\n")
	print(round(stats, 4))

	if(!is.na(x$method) & x$method != "chained" &
		x$type != "circle-arc")
		print(round(x$synthstats, 4))
	cat("\n")

	if(x$type != "equipercentile") {
		cat("Coefficients:\n")
		print(round(x$coef, 4))
		cat("\n")
	}

	invisible(x)
}
