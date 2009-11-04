\name{px}
\Rdversion{1.1}
\alias{px}
\title{Percentile Rank Function}
\description{
This function returns the percentile ranks for each score point in a frequency distribution.
}
\usage{
px(x)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point
}
}
\value{
Returns a vector of percentile ranks
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{fx}}
}
\examples{
x <- freqtab(ACTmath[,2],0:40)
px(x)
}
\keyword{distribution}