\name{freqtab}
\alias{freqtab}
\alias{as.freqtab}
\alias{is.freqtab}
\alias{mfreqtab}
\title{Frequency Distribution Tables}
\description{
Functions for creating, converting, and testing for frequency tables of class \dQuote{\code{freqtab}}.
}
\usage{
freqtab(x, v, xitems = 1:ncol(x), vitems = 1:ncol(v),
	xscale, vscale, na.rm = TRUE)

as.freqtab(x)

is.freqtab(x)

mfreqtab(x, margin = 1)
}
\arguments{
  \item{x, v}{
either vectors of total scores, one per examinee, or matrices of item responses with which total scores will be calculated. \code{x} refers to the total test, and \code{v} refers to the anchor test. See below for details.
}
  \item{xitems, vitems}{
vectors of column numbers or names with which total scores will be computed for \code{x} and \code{v} when they contain item responses.
}
  \item{xscale, vscale}{
the score scales for forms X and V, as vectors.
}
  \item{na.rm}{
logical with default \code{TRUE} specifying whether or not missing item responses should be ignored, i.e., treated as 0, when calculating total scores.
}
  \item{margin}{
integer specifying the margin over which frequencies should be summed when going from a bivariate to a univariate distribution.
}
}
\details{
\code{freqtab} creates a frequency table for a vector, or vectors, of scores. \code{as.freqtab} converts to class \dQuote{\code{freqtab}}, and \code{is.freqtab} tests for it. \code{mfreqtab} converts a frequency table from bivariate to univariate, but summing over margin \code{margin}.

With \code{freqtab}, if \code{x} and \code{v} are vectors of total scores, one per examinee, they will be factored with \code{levels = xscale}. This ensures that unobtained scores (i.e., scores not in \code{x} or \code{v}) will still be represented in the frequency tabulation.

When \code{x} is a matrix with at least two columns, total scores are calculated for each row across the columns specified in \code{xitems}. The same is true for \code{v} and \code{vitems}. If \code{x} and \code{vitems} are supplied but \code{v} is not, \code{x} is assumed to contain item responses for the anchor items as well.

\code{as.freqtab} serves to simply change a frequency table \code{x} to class \dQuote{\code{freqtab}}, as required by \code{\link{equate}} and other functions. The first column in \code{x} should contain the score scale. For a univariate distribution, the second column should contain the numbers of examinees obtaining each score in column 1. For a bivariate distribution, the second column in \code{x} should contain the anchor test score scale, and the third column the number of examinees obtaining each score combination.
}
\value{
In the univariate case (i.e., form V is not specified) a two-column matrix is returned, where column 1 contains the score scale and column 2 the frequency for each scale value. In the bivariate case a three-column matrix is returned where columns 1 and 2 contain all combinations of form X score scale points (in column 1) and score scale points on the anchor test (in column 2), and column 3 contains the total number of examinees obtaining the corresponding score combination.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{table}}, \code{\link{summary.freqtab}}
}
\examples{
set.seed(2005)
x <- round(rnorm(1000, 100, 10))
freqtab(x, xscale = 70:130)

as.freqtab(ACTmath[, 1:2])

y <- freqtab(x = KBneat$y[, 1], v = KBneat$y[, 2],
  xscale = 0:36, vscale = 0:12)
mfreqtab(y, margin = 2)
}
\keyword{methods}
