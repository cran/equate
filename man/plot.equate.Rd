\name{plot.equate}
\alias{plot.equate}
\title{Plotting of Equating Functions}
\description{
This function plots equating functions from one or more objects of class \dQuote{\code{equate}}.
}
\usage{
\method{plot}{equate}(..., add = FALSE, addident = TRUE, xpoints,
	ypoints, xlab = "Raw Scale", ylab = "Equated Scale",
	main = "Equating Plot", col = rainbow(length(x)),
	lty = 1, pchcol = "gray", addlegend = TRUE,
	legendtext, morepar = NULL)
}
\arguments{
  \item{\dots}{
One or more equating objects, each containing results for equating the same two test forms.
}
  \item{add}{
logical, with default \code{FALSE}, specifying whether to create a new plot or add to the current one.
}
  \item{addident}{
logical, with default \code{TRUE}, for plotting the identity function
}
  \item{xpoints, ypoints}{
optional vectors of the same length containing raw scores on forms X and Y, assuming an equivalent-groups design
}
  \item{xlab, ylab, main, col, lty}{
graphical parameters passed to \code{par}, with the lengths of col and lty recycled as necessary
}
  \item{pchcol}{
color used when plotting \code{xpoints} and \code{ypoints}
}
  \item{addlegend}{
logical, with default \code{TRUE}, indicating whether or not a legend should be added
}
  \item{legendtext}{
character vector of text to be passed to the \code{legend} argument of the \code{legend} function, defaulting to a combination of the equating types and methods specified in each equating object
}
  \item{morepar}{
additional graphical parameters passed to \code{par}
}
}
\details{
Equating functions are plotted as lines based on the concordance table for each equating object that is supplied. The resulting plot can be used to compare each function across the score scale.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\examples{
nx <- freqtab(KBneat$x[, 1], xscale = 0:36,
	v = KBneat$x[, 2], vscale = 0:12)
ny <- freqtab(KBneat$y[, 1], xscale = 0:36,
	v = KBneat$y[, 2], vscale = 0:12)

neq1 <- equate(nx, ny, type = "Linear", method = "Tuck")
neq2 <- equate(nx, ny, type = "equip", method = "freq")
neq3 <- equate(nx, ny, type = "linear", method = "braun")

plot(neq1, neq2, neq3)
}
\keyword{misc}
