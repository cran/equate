\name{descriptives}
\alias{summary.freqtab}
\alias{descript}
\alias{mean.freqtab}
\alias{sd.freqtab}
\alias{skew.freqtab}
\alias{kurt.freqtab}
\alias{max.freqtab}
\alias{min.freqtab}
\alias{cov.freqtab}
\alias{cor.freqtab}
\alias{px}
\alias{fx}
\title{Descriptive Statistics for Frequency Tables}
\description{
These functions return descriptive statistics for a frequency table of class \dQuote{\code{\link{freqtab}}}.
}
\usage{
\method{summary}{freqtab}(object, ...)

\method{mean}{freqtab}(x, ...)

sd.freqtab(x)

skew.freqtab(x)

kurt.freqtab(x)

\method{min}{freqtab}(x, ..., na.rm = FALSE)

\method{max}{freqtab}(x, ..., na.rm = FALSE)

cov.freqtab(x)

cor.freqtab(x)
}
\arguments{
  \item{object, x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point. If \code{x} is a three-column bivariate frequency table, \code{summary}, \code{cov.freqtab}, and \code{cor.freqtab} will utilize both variables, otherwise, only the first is used (see examples).
}
  \item{\dots}{
further arguments passed to or from other methods.
}
  \item{na.rm}{
logical indicating whether missing values should be removed, currently ignored since frequency tables cannot contain missing values.
}
}
\details{
\code{mean.freqtab}, \code{sd.freqtab}, \code{skew.freqtab}, and \code{kurt.freqtab} return the mean, standard deviation, skewness, and kurtosis. \code{min.freqtab} and \code{max.freqtab} return the maximum and minimum observed scores. \code{cov.freqtab} and \code{cor.freqtab} return the covariance and correlation for a bivariate frequency table.
}
\value{
\code{summary} returns a data frame of summary statistics, including the mean, standard deviation, skewness, kurtosis, minimum, maximum, and number of observations for each variable in \code{x}. Otherwise, a vector of length one is returned with the corresponding statistic.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}
}
\examples{
x <- as.freqtab(ACTmath[, 1:2])
summary(x)

y <- freqtab(KBneat$y[, 1], KBneat$y[, 2],
  xscale = 0:36, vscale = 0:12)
summary(y)
cov.freqtab(y)
cov(KBneat$y[, 1], KBneat$y[, 2])
}
\keyword{univar}
