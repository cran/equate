\name{se.eq}
\Rdversion{1.1}
\alias{se.eq}
\title{Standard Error of Equipercentile Equating}
\description{
This function returns the standard error of equipercentile equating, under the random groups design, using an equation derived by Lord (1982, p. 168).
}
\usage{
se.eq(q, g0, gm, xn, yn)
}
\arguments{
  \item{q, g0, gm, xn, yn}{
the five components required in the computation (for definitions see Lord, 1982)
}
}
\details{
This function is called by \code{\link{equate.eq}}, where it is only applicable to the random groups equating design
}
\value{
Returns the standard error for an equated score
}
\references{
Lord, F. (1982). The standard error of equipercentile equating. \emph{Journal of Educational Statistics, 7}(3), 165-174.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{equate.eq}}
}
\keyword{internal}