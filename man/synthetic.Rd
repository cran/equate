\name{synthetic}
\Rdversion{1.1}
\alias{synthetic}
\title{Synthetic Population Values}
\description{
This function returns the score distribution for a hypothetical synthetic examinee population, as required when equating under the nonequivalent groups design.
}
\usage{
synthetic(x, y, w, method, internal = TRUE)
}
\arguments{
  \item{x, y}{
bivariate score distributions of class \dQuote{\code{\link{freqtab}}} for forms X, with the anchor, and Y with the anchor. Columns 1 and 2 of each include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
}
  \item{w}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w})
}
  \item{method}{
string specifying the equating method, where the options are \code{"tucker"}, \code{levine}, and \code{"frequency"} (for frequency estimation). For the first two, equating \code{type=linear} is assumed, which can only be overridden via \code{\link{equate}} and \code{\link{equate.ln}}. Frequency estimation implies \code{type=equipercentile}
}
  \item{internal}{
logical indicating whether or not the anchor item scores are included in the total scores (default is \code{TRUE}). This is only required for the Levine method
}
}
\details{
The synthetic population score distribution is based on a weighted combination of the form X and Y distributions. The concept of a synthetic population was first described by Braun and Hollan (1982) and Kolen and Brennan (2004) summarize the assumptions made by each equating method.
}
\value{
For the Tucker and Levine methods, a list of length 1 (\code{$synthstats}) including the mean and standard deviation for form X, form Y, form XV, form YV, and for the synthetic popultion taking forms X and Y.

For the frequency estimation method, a frequency table for the synthetic population taking forms X and Y (\code{$synthtab}) is also returned.
}
\references{
Braun, H. I., & Holland, P. W. (1982). Observed-score test equating: A mathematical analysis of some ETS equating procedures. In P. W. Holland and D. B. Rubin (Eds.), \emph{Test Equating} (pp. 9-49). New York: Academic.

Kolen, M. J., & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\examples{
x <- KBneat$x
x <- freqtab(x[,1],0:36,x[,2],0:12)
y <- KBneat$y
y <- freqtab(y[,1],0:36,y[,2],0:12)
synthetic(x,y,w=1,method="T")
synthetic(x,y,w=1,method="L")
synthetic(x,y,w=1,method="F")
}
\keyword{methods}