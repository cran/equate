\name{freqavg}
\Rdversion{1.1}
\alias{freqavg}
\title{Frequency Averaging}
\description{
This function replaces frequencies falling below a specified threshold with an average of adjacent values. Frequency averaging is a simple method for smoothing a score distribution, typically one with unobserved scores (for a description see Moses & Holland, 2008).
}
\usage{
freqavg(x, jmin = 1)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point
}
  \item{jmin}{
the minimum frequency, as an integer, below which frequencies will be replaced
}
}
\value{
Returns a vector of adjusted frequencies. If no values fall below \code{jmin} these will simply be the original frequencies.
}
\references{
Moses, T., & Holland, P. (2008). \emph{Notes on a general framework for observed score equating} (ETS Research Rep. No. RR-08-59). Princeton, NJ: ETS.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqbump}}, \code{\link{loglinear}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)

# compare freqavg to freqbump:
xavg <- freqavg(xtab)
xbump <- freqbump(xtab[,2],Kx=60)
cbind(xtab,avg=xavg,bump=round(xbump*sum(xtab[,2]),6))
}
\keyword{smooth}
