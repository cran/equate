\name{freqadjust}
\alias{freqavg}
\alias{freqbump}
\title{Frequency Adjusting and Averaging}
\description{
These functions provide simple methods for smoothing score distributions, typically ones with unobserved score points, prior to performing equipercentile equating. \code{freqbump} smooths by adding a small relative frequency to each score while adjusting the probabilities to sum to one (see Kolen and Brennan, 2004, p. 48). \code{freqavg} replaces frequencies falling below a specified threshold with an average of adjacent values (see Moses and Holland, 2008).
}
\usage{
freqbump(x, jmin = 1e-6, Kx = max(x[,1]))
freqavg(x, jmin = 1)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point. For \code{freqbump}, \code{x} can also be a bivariate frequency table.
}
  \item{jmin}{
for \code{freqbump}, the frequency to be added to each score point (default is 1e-6)
for \code{freqavg}, the minimum frequency, as an integer, below which frequencies will be replaced (default is 1)
}
  \item{Kx}{
integer indicating the number of items in the scale
}
}
\value{
Returns a vector of adjusted frequencies. For \code{freqavg}, if no values fall below \code{jmin} these will simply be the original frequencies. For \code{freqbump}, frequencies are in the form of probabilities, summing to 1. 
}
\references{
Kolen, M. J., and Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.

Moses, T., and Holland, P. W. (2008). \emph{Notes on a general framework for observed score equating} (ETS Research Rep. No. RR-08-59). Princeton, NJ: ETS.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{loglinear}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000, 100, 10))
xtab <- freqtab(x, xscale = 70:130)

xavg <- freqavg(xtab)
xbump <- freqbump(xtab)
cbind(xtab, avg = xavg,
  bump = round(xbump * sum(xtab[, 2]), 6))
}
\keyword{smooth}
