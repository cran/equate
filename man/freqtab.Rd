\name{freqtab}
\Rdversion{1.1}
\alias{freqtab}
\title{Frequency Tabulation}
\description{
This function creates a frequency table for a vector, or vectors, of scores.
}
\usage{
freqtab(x, xscale, v, vscale, addclass = FALSE)
}
\arguments{
  \item{x}{
vector of scores (one per examinee) for form X
}
  \item{xscale}{
the score scale for form X, as a vector
}
  \item{v}{
vector of scores (one per examinee) for form X examinees on the anchor test
}
  \item{vscale}{
the anchor test score scale, as a vector
}
  \item{addclass}{
logical. When \code{TRUE} (and \code{v} is missing) \code{x} and \code{xscale} are simply combined as a matrix of class \dQuote{\code{freqtab}}
}
}
\details{
If \code{x} is a vector of total scores, values are factored with \code{levels = xscale}. This ensures that unobtained scores (i.e., scores not in \code{x}) will still be represented in the tabulation. If \code{x} is a vector of frequencies, with one value for each score scale point provided in \code{xscale}, and \code{addclass=TRUE}, the function serves to simply change the class of \code{x} to \dQuote{\code{freqtab}}.
}
\value{
In the univariate case (i.e., no anchor is specified) a 2-column matrix is returned, where column 1 contains the score scale and column 2 the frequency for each scale value. In the bivariate case a 3-column matrix is returned where columns 1 and 2 contain all combinations of form X score scale points (in column 1) and score scale points on the anchor test (in column 2), and column 3 contains the total number of examinees obtaining the corresponding score combination.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{table}}
}
\examples{
set.seed(2005)
x1 <- round(rnorm(1000,100,10))
x1scale <- 70:130
freqtab(x1,x1scale)

x2 <- ACTmath[,2]
x2scale <- ACTmath[,1]
freqtab(x2,x2scale,addclass=TRUE)

xscale <- 0:36
vscale <- 0:12
freqtab(KBneat$x[,1],xscale,KBneat$x[,2],vscale)
}
\keyword{distribution}
