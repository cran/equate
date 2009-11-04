\name{kurt.freqtab}
\Rdversion{1.1}
\alias{kurt.freqtab}
\title{Kurtosis of a Frequency Distribution}
\description{
This function returns the kurtosis of a frequency table of class \dQuote{\code{\link{freqtab}}}
}
\usage{
kurt.freqtab(x)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point
}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{skew.freqtab}}, \code{\link{cov.freqtab}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)
kurt.freqtab(xtab)
}
\keyword{univar}
