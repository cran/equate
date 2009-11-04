\name{se.boot}
\Rdversion{1.1}
\alias{se.boot}
\title{Bootstrap Standard Errors of Equating}
\description{
This function returns the standard deviation of equated scores over multiple replications of a specified equating.
}
\usage{
se.boot(x, y, xn = sum(x[,ncol(x)]), yn = sum(y[,ncol(y)]), 
  reps = 500, eqfun, ...)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y. Under the random groups design (i.e., no anchor test) \code{x} and \code{y} will each contain the score scale in column 1, and the number of examinees obtaining each score in column 2. For the nonequivalent groups design bivariate frequency tables are used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
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