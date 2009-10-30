\name{skew.freqtab}
\Rdversion{1.1}
\alias{skew.freqtab}
\title{Skewness of a Frequency Distribution}
\description{
This function returns the skewness of a frequency table of class \dQuote{\code{\link{freqtab}}}
}
\usage{
skew.freqtab(freqtab)
}
\arguments{
  \item{freqtab}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point
}
}
\value{
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{kurt.freqtab}}, \code{\link{var.freqtab}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)
skew.freqtab(xtab)
}
\keyword{univar}