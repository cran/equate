\name{plot.freqtab}
\alias{plot.freqtab}
\title{Plotting of Frequency Distributions}
\description{
This function plots univariate and bivariate frequency tables of class \dQuote{\code{\link{freqtab}}}
}
\usage{
\method{plot}{freqtab}(x, y, col1 = "gray27", col2 = "gray74",
	type = "h", pch=16, yoffset = 0.3, ...)
}
\arguments{
  \item{x}{
univariate or bivariate score distribution of class \dQuote{\code{\link{freqtab}}}.
}
  \item{y}{
univariate score distribution, either as an object of class \dQuote{\code{\link{freqtab}}} or a vector of frequencies, to be added to the plot of univariate \code{x} (see details).
}
  \item{col1}{
color of lines used in plotting univarite \code{x}, and in plotting both variables in bivariate \code{x}.
}
  \item{col2}{
color of lines used in plotting \code{y}.
}
  \item{type}{
character vector describing the type of plot to be drawn when \code{x} is univariate, where \code{type[1]} is used for plotting \code{x} and \code{type[2]} for plotting \code{y}. The default \code{"h"} plots a vertical line for each frequency. \code{type} is recycled if necessary.
}
  \item{pch}{
plotting point (see \code{\link{par}}) used in bivariate scatterplot.
}
  \item{yoffset}{
integer indicating amount of offset when plotting vertical lines (i.e., \code{type[2] == "h"}) for \code{y}, expressed in terms of \code{x} scale points.
}
  \item{\dots}{
further arguments passed to or from other methods, such as graphical parameters besides \code{col}.
}
}
\details{
When \code{x} is bivariate, e.g., coming from one form of a nonequivalent groups equating design, \code{y} is ignored. Otherwise \code{x} and \code{y} are assumed to be on similar scales, where \code{x} determines the x-axis.

Colors must be specified using \code{col1} and \code{col2}.
}
\value{
The univariate option produces a single, simple, line plot of \code{type="h"}. The bivariate option produces a scatterplot with a marginal barplot for each distribution.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{plot.table}}, \code{\link{lines}}
}
\examples{
neat.xv <- freqtab(KBneat$x[, 1], KBneat$x[, 2],
	xscale = 0:36, vscale = 0:12)
plot(neat.xv)

neat.x <- freqtab(KBneat$x[, 1], xscale = 0:36)
neat.v <- freqtab(KBneat$x[, 2], vscale = 0:12)
plot(neat.x, neat.v, lwd = 2)
plot(neat.x, neat.v, lwd = 2, type = "l")
}
\keyword{misc}
