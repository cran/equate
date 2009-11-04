\name{mean.freqtab}
\Rdversion{1.1}
\alias{mean.freqtab}
\title{Mean of a Frequency Table}
\description{
This function returns the mean of a frequency table of class \dQuote{\code{\link{freqtab}}}
}
\usage{
mean.freqtab(x, ...)
}
\arguments{
  \item{x}{
two-column matrix of class \code{freqtab}, where column 1 specifies the score scale and column 2 the frequencies for each score point. If \code{x} is a bivariate frequency table, only the first variable will be considered
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}, \code{\link{cov.freqtab}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000,100,10))
xscale <- 70:130
xtab <- freqtab(x,xscale)
mean.freqtab(xtab); mean(x)

yscale <- 0:36
vscale <- 0:12
ny <- freqtab(KBneat$y[,1],yscale,KBneat$y[,2],vscale)
mean.freqtab(ny); mean(KBneat$y[,1])
}
\keyword{univar}