\name{equate.ln}
\Rdversion{1.1}
\alias{equate.ln}
\title{Linear Equating}
\description{
This function (called by \code{\link{equate}}) conducts linear equating of test scores under the random groups and nonequivalent groups with anchor test designs. Types of equating include mean and linear, and methods include Tucker and Levine.
}
\usage{
equate.ln(x, y, type = "linear", method = "none", w = 1, 
  internal = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y. Under the random groups design (i.e., no anchor test) \code{x} and \code{y} will each contain the score scale in column 1, and the number of examinees obtaining each score in column 2. For the nonequivalent groups design bivariate frequency tables are used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
}
  \item{type}{
the type of equating, implying the numer of coefficients used in the linear transformation, where \code{"mean"} relies on the intercept and \code{"linear"} (default) the intercept and slope
}
  \item{method}{
string specifying the equating method, where \code{"none"} (default) indicates observed score equating under the random groups design, and both \code{"tucker"} and \code{"levine"} indicate the respective methods under the nonequivalent groups design
}
  \item{w}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w}) when estimating the synthetic population
}
  \item{internal}{
logical indicating whether or not the anchor item scores are included in the total scores (default is \code{TRUE})
}
  \item{verbose}{
logical indicating whether output beyond a vector of equated scores should be supplied (\code{TRUE}) or not (\code{FALSE}, default)
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
For the random groups design \code{x} and \code{y} are separate frequency tables where each row contains 1) a score scale point and 2) the total number of examinees obtaining that score. For the nonequivalent groups design each row contains 1) a score scale point on \code{x} or \code{y}, 2) a score scale point on the anchor test, and 3) the total number of examinees obtaining that score combination (see \code{\link{freqtab}} for more details).

Two equating types (mean and linear) and two equating methods (Tucker and Levine) are currently supported. Varying the \code{type} and \code{method} arguments results in four possible equating scenarios. Whenever \code{"tucker"} or \code{"levine"} is specified, the nonequivalent groups with anchor test design is assumed. In this case anchor test scores must be included in \code{x} and \code{y}, so that each contains 3 columns, and the argument \code{w} is also required (see \code{\link{synthetic}} for details).
}
\value{
When \code{verbose = FALSE} a vector of equated scores is returned, otherwise a list is returned containing the following components:
  \item{yx}{form Y equivalents of the form X scores. For the random groups design, when \code{verbose = TRUE} a column of standard errors is included}
  \item{coefficients}{the slope and intercept used in the linear conversion}
  \item{synthstats}{means and standard deviations for the anchor and synthetic distributions}
  \item{anchortab}{frequency table for the anchor test}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{\code{\link{equate}}, \code{\link{synthetic}}}
\examples{
# see the function 'equate' for more examples
x <- freqtab(ACTmath[,2],0:40,addclass=TRUE)
y <- freqtab(ACTmath[,3],0:40,addclass=TRUE)
cbind(0:40,equate.ln(x,y,type="m"),
  equate.ln(x,y,type="l"))
}
\keyword{methods}