\name{synthetic}
\Rdversion{1.1}
\alias{synthetic}
\title{Synthetic Population Values}
\description{
This function returns the score distribution for a hypothetical synthetic examinee population, as required when equating under the nonequivalent groups design.
}
\usage{
synthetic(x, xv, y, yv, w, method, internal = TRUE, scale, vscale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of total scores for form X (one per examinee)
}
  \item{xv}{
vector of scores (one per examinee) for form X examinees on the anchor test
}
  \item{y}{
vector of total scores for form Y (one per examinee)
}
  \item{yv}{
vector of scores (one per examinee) for form Y examinees on the anchor test
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
  \item{scale}{
the score scale, as a vector (only required for frequency estimation)
}
  \item{vscale}{
the score scale for the anchor test, as a vector (only required for frequency estimation)
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
x <- KBneat$x[,1]
xv <- KBneat$x[,2]
y <- KBneat$y[,1]
yv <- KBneat$y[,2]
synthetic(x,xv,y,yv,w=1,method="T")
synthetic(x,xv,y,yv,w=1,method="L")
synthetic(x,xv,y,yv,w=1,method="F",scale=0:36,vscale=0:12)
}
\keyword{multivariate}