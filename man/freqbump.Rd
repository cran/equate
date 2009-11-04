\name{freqbump}
\Rdversion{1.1}
\alias{freqbump}
\title{Frequency Adjustment}
\description{
This function smooths a frequency distribution by adding a small relative frequency to each score while adjusting the probabilities to sum to one (see Kolen & Brennan, 2004, p. 48).
}
\usage{
freqbump(x, jmin = 10^-6, Kx = length(x))
}
\arguments{
  \item{x}{
vector of either counts (summing to the number of examinees) or probabilities (summing to 1). If \code{x} does not sum to one it is assumed to be a vector of counts and it is converted to probabilities
}
  \item{jmin}{
the frequency to be added to each score point (default is 10^-6)
}
  \item{Kx}{
integer indicating the number of items in the scale (default is \code{length(x)})
}
}
\value{
Returns a vector of adjusted relative frequencies, as probabilities
}
\references{
Kolen, M. J., & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}, \code{\link{freqavg}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)

# compare freqbump to freqavg:
xbump <- freqbump(xtab[,2])
xavg <- freqavg(xtab)
cbind(xtab,bump=round(xbump*sum(xtab[,2]),6),avg=xavg)
}
\keyword{smooth}
