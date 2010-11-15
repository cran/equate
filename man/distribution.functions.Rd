\name{percentiles}
\alias{px}
\alias{fx}
\title{Percentile Rank and Cumulative Distribution Functions}
\description{
\code{px} returns the percentile ranks and \code{fx} the cumulative frequencies for each score point in a frequency distribution.
}
\usage{
px(x,y)
fx(x)
}
\arguments{
  \item{x}{
two-column matrix of class \dQuote{\code{\link{freqtab}}}, where column 1 specifies the score scale and column 2 the frequencies for each score point. For \code{px}, \code{x} can also be a vector of score points (see details)
}
  \item{y}{
optional frequency distribution of class \dQuote{\code{\link{freqtab}}} (see details)
}
}
\details{
\code{px} can be used to find the percentile ranks for a distribution of scores. When \code{x} is a frequency table and \code{y} is missing, the X distribution is used. If a frequency table is supplied for \code{y}, the percentile ranks are found for each score from \code{x} in the Y distribution. In this latter case \code{x} can either be a frequency table or simply a vector of scores.
}
\value{
\code{px} returns a vector of percentile ranks, and \code{fx} a vector of cumulative frequencies.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\examples{
x <- freqtab(0:40, ACTmath[, 2])
px(x)
fx(x)
px(c(8.2, 40, 41), x)
}
\keyword{distribution}
