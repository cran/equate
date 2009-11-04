\name{cov.freqtab}
\Rdversion{1.1}
\alias{cov.freqtab}
\title{Covariance of a Frequency Table}
\description{
This function computes the covariance of a frequency table of class \dQuote{\code{\link{freqtab}}}
}
\usage{
cov.freqtab(x)
}
\arguments{
  \item{x}{
score distribution of class \dQuote{\code{\link{freqtab}}} for form X. Under the random groups design (i.e., no anchor test) \code{x} will contain the score scale in column 1, and the number of examinees obtaining each score in column 2. For the nonequivalent groups design a bivariate frequency table is used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
}
}
\value{
Returns the variance, in the univariate case, and the covariance for the bivariate
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
cov.freqtab(xtab); var(x)

yscale <- 0:36
vscale <- 0:12
ny <- freqtab(KBneat$y[,1],yscale,KBneat$y[,2],vscale)
cov.freqtab(ny); cov(KBneat$y[,1],KBneat$y[,2])
}
\keyword{univar}