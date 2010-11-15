\name{plot.freqtab}
\alias{plot.freqtab}
\title{Plotting of Frequency Distributions}
\description{
This function plots univariate and bivariate frequency tables of class \dQuote{\code{\link{freqtab}}}
}
\usage{
plot.freqtab(x, y, col1 = "gray27", col2 = "gray74", pch=16,
  yoffset = 0.3, ...)
}
\arguments{
  \item{x}{
univariate or bivariate score distribution of class \dQuote{\code{\link{freqtab}}}
}
  \item{y}{
univariate score distribution, also of class \dQuote{\code{\link{freqtab}}}, to be added to the plot of univariate \code{x} (see details)
}
  \item{col1}{
color of lines used in plotting univarite \code{x}, and in plotting both variables in bivariate \code{x}
}
  \item{col2}{
color of lines used in plotting \code{y}
}
  \item{pch}{
plotting point (see \code{\link{par}}) used in bivariate scatterplot
}
  \item{yoffset}{
integer indicating amount of offset, in terms of \code{x} scale points, when plotting lines for \code{y}
}
  \item{\dots}{
further arguments passed to or from other methods, such as graphical parameters besides \code{col}
}
}
\details{
When \code{x} is bivariate \code{y} is ignored. Otherwise \code{x} and \code{y} are (currently) assumed to be on the same scale, which determines the x-axis.

Colors must be specified using \code{col1} and \code{col2}
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
xscale <- 0:36
vscale <- 0:12
neat.xv <- freqtab(xscale, KBneat$x[, 1],
  vscale, KBneat$x[, 2])
plot(neat.xv)

neat.x <- freqtab(xscale, KBneat$x[, 1])
neat.v <- freqtab(vscale, KBneat$x[, 2])
plot(neat.x, neat.v, lwd = 2)
}
\keyword{misc}
