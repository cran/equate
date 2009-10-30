\name{equate.eq}
\Rdversion{1.1}
\alias{equate.eq}
\title{Equipercentile Equating}
\description{
This function (called by \code{\link{equate}}) conducts equipercentile equating of test scores under the random groups and nonequivalent groups with anchor test designs (see Kolen & Brennan, 2004). Methods for presmoothing scores are supported.
}
\usage{
equate.eq(x, y, scale, method = "none", Ky = max(scale), 
  xv, yv, vscale, w, smooth = "none",  jmin, xscorefun, 
  yscorefun, verbose = FALSE, ...)
}
\arguments{
  \item{x, y}{
vectors of scores for forms X and Y, either as total scores (one per examinee), or total counts (one per score scale point). Total scores are required for the nonequivalent groups design
}
  \item{scale}{
the score scale as a vector
}
  \item{method}{
string specifying the equating method, where \code{"none"} (default) indicates observed score equating under the random groups design, and \code{"frequency"} indicates frequency estimation with nonequivalent groups and requires additional anchor test scores
}
  \item{Ky}{
integer, typically indicating the number of items on form Y, or the maximum possible score on Y. This value is used in equating the extremes of the scale (see below for more details)
}
  \item{xv, yv}{
vectors of anchor test total scores for form X and Y examinees
}
  \item{vscale}{
the score scale for the anchor test, as a vector
}
  \item{w}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w})
}
  \item{smooth}{
string indicating one of four smoothing methods: \code{"none"} (default), \code{"bump"}, \code{"average"}, and \code{"loglin"} (see below for details)
}
  \item{jmin}{
smoothing parameter required for methods \code{"bump"} and \code{"average"}
}
  \item{xscorefun, yscorefun}{
matrices, required by smoothing method \code{"loglin"}, which specify the score functions for forms X and Y (see below for details)
}
  \item{verbose}{
logical indicating whether output beyond a vector of equated scores should be supplied (\code{TRUE}) or not (\code{FALSE}, default)
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
As noted above, for the random groups design \code{x} and \code{y} can be supplied as either total scores, where the length of each will equal the number of examinees, or as counts, where the length of each will equal the length of the score scale, with each value indicating the total number of examinees obtaining the corresponding score. With the nonequivalent groups design, only total scores may be used.

Two equating methods are currently supported: observed score equipercentile equating, which assumes a random groups equating design, and frequency estimation equipercentile equating, which assumes a nonequivalent groups with anchor test design. Smoothing via \code{equate.eq} is only possible for the random groups design, where the additional argument \code{jmin} is required for methods \code{"bump"} and \code{"average"} (see \code{\link{freqbump}} and \code{\link{freqavg}}). Loglinear smoothing requires additional parameters (see \code{\link{loglinear}}). If smoothed distributions are desired for the frequency estimation method they must be supplied. The frequency estimation method also requires the additional arguments \code{xv}, \code{yv}, \code{vscale}, and \code{w}, which are passed to \code{\link{synthetic}}.

The parameter \code{Ky} is used to obtain form Y equivalents of form X scores with percentile ranks of 100. Typically \code{Ky} is set to be the number of score points in the form Y scale, which assumes that scores are integers ranging from 1 (or 0) to the total number of items, and that each item is scored correct/incorrect. Scores on other scales (such as scales which include negative values, or which exclude zero) may also be used. In such cases \code{Ky} be set to the highest possible score on form Y, or alternatively the highest observed score on Y. For optimal results, scales should be monotonically increasing integer sequences. As another alternative, \code{Ky = "i"} will leave all form X scores with percentiles of 0 and 100 in raw form (i.e., unequated).
}
\value{
Returns a list including the following components when \code{verbose = TRUE} (only the first is returned when \code{verbose = FALSE}):
  \item{yx}{form Y equivalents of the form X scores. For the random groups design, and when \code{verbose = TRUE}, a column of standard errors is included}
  \item{smoothmethod}{smoothing method}
  \item{smoothout}{output corresponding to the particular smoothing method (see smoothing functions for details)}
  \item{synthstats}{matrix of means and standard deviations for the anchor and synthetic distributions}
  \item{synthtab}{frequency table for the synthetic distribution}
  \item{anchortab}{frequency table for the anchor test}
}
\references{
Kolen, M. J., & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\note{}
\seealso{\code{\link{equate}}, \code{\link{synthetic}}, \code{\link{freqbump}}, \code{\link{freqavg}}, \code{\link{loglinear}}}
\examples{
# see the function 'equate' for more examples
x <- KBneat$x
y <- KBneat$y
scale <- 0:36
equate.eq(x[,1],y[,1],scale,method="f",
  xv=x[,2],yv=y[,2],w=.5,vscale=0:12)
}
\keyword{methods}