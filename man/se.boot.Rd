\name{se.boot}
\Rdversion{1.1}
\alias{se.boot}
\title{Bootstrap Standard Errors of Equating}
\description{
This function calculates the standard deviation of equated scores over multiple replications of a specified equating.
}
\usage{
se.boot(x, y, xn = sum(x[, 2]), yn = sum(y[, 2]), reps = 500, eqfun, ...)
}
\arguments{
  \item{x, y}{
two-column matrices of class \code{freqtab}, where column 1 specifies the score scale and column 2 the frequencies for each score point
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
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
Samples are drawn of size \code{xn} and \code{yn}, with replacement, from each score distribution. Form Y equivalents of each form X score are then obtained using the specified equating function \code{eqfun}. This process is repeated \code{reps} times.

This function is called by \code{\link{equate}}.
}
\value{
Returns a vector of standard deviations, one per score scale point.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\examples{
xtab <- freqtab(KBneat$x[,1],0:36)
ytab <- freqtab(KBneat$y[,1],0:36)
se.boot(xtab,ytab,reps=100,eqfun="equate.ln")
}
\keyword{univar}