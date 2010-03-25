\name{descriptives}
\alias{descript}
\alias{mean.freqtab}
\alias{sd.freqtab}
\alias{skew.freqtab}
\alias{kurt.freqtab}
\alias{max.freqtab}
\alias{min.freqtab}
\alias{cov.freqtab}
\title{Descriptive Statistics for a Frequency Table}
\description{
The functions \code{mean.freqtab}, \code{sd.freqtab}, \code{skew.freqtab}, and \code{kurt.freqtab} return the mean, standard deviation, skewness, and kurtosis for a frequency table of class \dQuote{\code{\link{freqtab}}}. \code{descript} returns all four of these values along with the number of observations.

\code{min.freqtab} and \code{max.freqtab} return the minimum and maximum.

\code{cov.freqtab} returns the variance of a univariate frequency table and the covariance of a bivariate.
}
\usage{
mean.freqtab(x, ...)
sd.freqtab(x)
skew.freqtab(x)
kurt.freqtab(x)
descript(x)
min.freqtab(x, ..., na.rm = FALSE)
max.freqtab(x, ..., na.rm = FALSE)
cov.freqtab(x)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point. If \code{x} is a three-column bivariate frequency table, only the first distribution is used, except with \code{cov.freqtab} (see examples)
}
  \item{\dots}{
further arguments passed to or from other methods
}
  \item{na.rm}{
logical indicating whether or not missing values should be removed
}
}
\value{
\code{descript} returns a vector which includes the mean, standard deviation, skewness, kurtosis, and number of observations. Otherwise a vector of length one is returned.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}
}
\examples{
xtab <- freqtab(ACTmath[,1], ACTmath[,2])
descript(xtab)
mean(xtab)

set.seed(2008)
y <- round(rnorm(100, 100, 10))
ytab <- freqtab(70:130, y)
min(y); min(ytab)
max(y); max(ytab)
cov.freqtab(ytab); var(y)

yscale <- 0:36
vscale <- 0:12
ny <- freqtab(yscale, KBneat$y[,1],
  vscale, KBneat$y[,2])
cov.freqtab(ny); cov(KBneat$y[,1], KBneat$y[,2])
}
\keyword{univar}
