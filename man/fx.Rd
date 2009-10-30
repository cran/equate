\name{fx}
\Rdversion{1.1}
\alias{fx}
\title{Cumulative Distribution Function}
\description{
This function returns cumulative frequencies for each score point in a frequency distribution.
}
\usage{
fx(freqtab)
}
\arguments{
  \item{freqtab}{
two-column matrix of class \code{\link{freqtab}}, where column 1 specifies the score scale and column 2 the frequencies for each score point
}
}
\details{
}
\value{
Returns a vector of cumulative frequencies.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{px}}
}
\examples{
x <- freqtab(ACTmath[,2],0:40)
fx(x)
}
\keyword{distribution}
