\name{se.boot}
\alias{se.boot}
\title{Bootstrap Standard Errors of Equating}
\description{
This function returns the standard deviation of equated scores over multiple replications of a specified equating.
}
\usage{
se.boot(x, y, xn = sum(x[,ncol(x)]), yn = sum(y[,ncol(y)]), 
  reps = 100, eqfun, returnboots = FALSE, ...)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y. These may be univariate or bivariate distributions, as required by \code{\link{equate}}
}
  \item{xn, yn}{
integers specifying the number of scores to sample from each distribution at each replication (default is the total number of scores in each)
}
  \item{reps}{
number of replications
}
  \item{eqfun}{
string indicating the equating function
}
  \item{returnboots}{
logical (default is \code{FALSE}) indicating whether or not the matrix of equated scores should be returned
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
Samples are drawn of size \code{xn} and \code{yn}, with replacement, from each score distribution. Form Y equivalents of each form X score are then obtained using the specified equating function \code{eqfun}. This process is repeated \code{reps} times.

This function is most easily called from \code{\link{equate}}.
}
\value{
Returns either a vector of standard deviations, one per score scale point, or a matrix of equated scores with number of columns equal to \code{reps}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\examples{
xtab <- freqtab(0:36, KBneat$x[, 1])
ytab <- freqtab(0:36, KBneat$y[, 1])
se.boot(xtab, ytab, reps = 100, eqfun = "equate.ln")
}
\keyword{methods}
