\name{freqtab}
\alias{freqtab}
\alias{as.freqtab}
\title{Frequency Tabulation}
\description{
\code{freqtab} creates a frequency table for a vector, or vectors, of scores. \code{as.freqtab} converts to class \dQuote{\code{freqtab}}
}
\usage{
freqtab(xscale, x, vscale, v)
as.freqtab(xscale, x)
}
\arguments{
  \item{xscale}{
the score scale for form X, as a vector
}
  \item{x}{
vector of scores (one per examinee) on items in form X (see details for other options)
}
  \item{vscale}{
the score scale for form V, as a vector
}
  \item{v}{
vector of scores (one per examinee) for form X examinees on items in form V (i.e., the anchor test)
}
}
\details{
If \code{x} is a vector of scores, one per examinee, values are factored with \code{levels = xscale}. This ensures that unobtained scores (i.e., scores not in \code{x}) will still be represented in the tabulation. If \code{x} is a two-column matrix, \code{as.freqtab} serves to simply change its class to \dQuote{\code{freqtab}}, as required by \code{\link{equate}}. In this case the first column should contain the score scale and the second column the number of examinees obtaining each score. Alternatively, both \code{xscale} and \code{x} can be included in \code{as.freqtab}, with the same result.
}
\value{
In the univariate case (i.e., form V is not specified) a two-column matrix is returned, where column 1 contains the score scale and column 2 the frequency for each scale value. In the bivariate case a three-column matrix is returned where columns 1 and 2 contain all combinations of form X score scale points (in column 1) and score scale points on the anchor test (in column 2), and column 3 contains the total number of examinees obtaining the corresponding score combination.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{table}}
}
\examples{
set.seed(2005)
x1 <- round(rnorm(1000, 100, 10))
x1scale <- 70:130
freqtab(x1scale, x1)

x2 <- ACTmath[, 2]
x2scale <- ACTmath[, 1]
as.freqtab(x2scale, x2)

xscale <- 0:36
vscale <- 0:12
freqtab(xscale, KBneat$x[, 1],
  vscale, KBneat$x[, 2])
}
\keyword{utilities}
