\name{equate.eq}
\Rdversion{1.1}
\alias{equate.eq}
\title{Equipercentile Equating}
\description{
This function (called by \code{\link{equate}}) conducts equipercentile equating of test scores under the random groups and nonequivalent groups with anchor test designs. Methods for presmoothing scores are supported.
}
\usage{
equate.eq(x, y, method = "none", Ky = max(y[,1]), w = 1, smooth = "none", 
  jmin, xscorefun, yscorefun, verbose = FALSE, ...)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y. Under the random groups design (i.e., no anchor test) \code{x} and \code{y} will each contain the score scale in column 1, and the number of examinees obtaining each score in column 2. For the nonequivalent groups design bivariate frequency tables are used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
}
  \item{method}{
string specifying the equating method, where \code{"none"} (default) indicates observed score equating under the random groups design, and \code{"frequency"} indicates frequency estimation with nonequivalent groups and requires additional anchor test scores
}
  \item{Ky}{
integer, typically indicating the number of items on form Y, or the maximum possible score on Y. This value is used in equating the extremes of the scale (see below for more details)
}
  \item{w}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w}) when estimating the synthetic population
}
  \item{smooth}{
string indicating one of four smoothing methods: \code{"none"} (default), \code{"bump"}, \code{"average"}, and \code{"loglin"} (see below for details)
}
  \item{jmin}{
smoothing parameter required for methods \code{"bump"} and \code{"average"}
}
  \item{xscorefun, yscorefun}{
matrices, required by smoothing method \code{"loglin"} which specify the score functions for forms X and Y (see below for details)
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

Two equating methods are currently supported: observed score equipercentile equating, which assumes a random groups equating design, and frequency estimation equipercentile equating, which assumes a nonequivalent groups with anchor test design. For the random groups design, all three smoothing options are available: \code{"bump"} and \code{"average"} (see \code{\link{freqbump}} and \code{\link{freqavg}}) require the additional argument \code{jmin} and loglinear smoothing (\code{"loglin"}) requires the score functions \code{xscorefun} and \code{yscorefun} (see \code{\link{loglinear}}). For frequency estimation and the nonequivalent groups design, only smoothing methods \code{"freqbump"} and \code{"loglin"} are supported.

For frequency estimation a \dQuote{synthetic} distribution of scores is used to equate X to Y. These procedures are carried out by \code{\link{synthetic}}.

The parameter \code{Ky} is used to obtain form Y equivalents of form X scores with percentile ranks of 100. Typically \code{Ky} is set to be the number of score points in the form Y scale, which assumes that scores are integers ranging from 1 (or 0) to the total number of items, and that each item is scored correct/incorrect. Scores on other scales (such as scales which include negative values, or which exclude zero) may also be used. In such cases \code{Ky} can be set to the highest possible score on form Y, or alternatively the highest observed score on Y. For optimal results, scales should be monotonically increasing integer sequences. Specifying \code{Ky = "i"} will leave all form X scores with percentiles of 0 and 100 in raw form (i.e., unequated).
}
\value{
When \code{verbose = FALSE} a vector of equated scores is returned, otherwise a list is returned containing the following components:
  \item{yx}{form Y equivalents of the form X scores. For the random groups design, and when \code{verbose = TRUE}, a column of standard errors is included}
  \item{smoothmethod}{smoothing method}
  \item{smoothout}{output corresponding to the particular smoothing method (see smoothing functions for details)}
  \item{synthstats}{means and standard deviations for the anchor and synthetic distributions}
  \item{synthtab}{frequency table for the synthetic distribution}
  \item{anchortab}{frequency table for the anchor test}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{\code{\link{equate}}, \code{\link{synthetic}}, \code{\link{freqbump}}, \code{\link{freqavg}}, \code{\link{loglinear}}}
\examples{
# see the function 'equate' for more examples
x <- KBneat$x
x <- freqtab(x[,1],0:36,x[,2],0:12)
y <- KBneat$y
y <- freqtab(y[,1],0:36,y[,2],0:12)
equate.eq(x,y,method="f",w=.5,smooth="b")
}
\keyword{methods}