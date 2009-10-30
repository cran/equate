\name{var.freqtab}
\Rdversion{1.1}
\alias{var.freqtab}
\title{Variance of a Frequency Table}
\description{
This function computes the variance of a frequency table of class \dQuote{\code{\link{freqtab}}}
}
\usage{
var.freqtab(freqtab)
}
\arguments{
  \item{freqtab}{
two-column matrix of class \dQuote{\code{\link{freqtab}}} where column 1 specifies the score scale and column 2 the frequencies for each score point
}
}
\details{
}
\value{
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}, \code{\link{mean.freqtab}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)
var.freqtab(xtab); var(x)
}
\keyword{univar}
