print.equate <- function(x, ...) {
  up <- function(x)
    gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  design <- up(x$design)
  if(x$method == "none")
    method <- up(x$type)
  else method <- paste(up(x$method), up(x$type))

  cat("\n", method, " Equating: ", design, sep = "")
  cat("\n\nSummary Statistics:\n")
  print(round(x$stats, 4))
  if(x$method != "none" & x$method != "chained" &
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
