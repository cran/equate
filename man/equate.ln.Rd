\name{equate.ln}
\Rdversion{1.1}
\alias{equate.ln}
\title{Linear Equating}
\description{
This function (called by \code{\link{equate}}) conducts linear equating of test scores under the random groups and nonequivalent groups with anchor test designs. Types of equating include mean and linear, and methods include Tucker and Levine (see Kolen & Brennan, 2004).
}
\usage{
equate.ln(x, y, scale, type = "linear", method = "none", 
  xv, yv, vscale, w = 1, internal = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{x, y}{
vectors of scores for forms X and Y, either as total scores (one per examinee), or total counts (one per score scale point). Total scores are required for the nonequivalent groups design
}
  \item{scale}{
the score scale, as a vector
}
  \item{type}{
the type of equating, implying the numer of coefficients used in the linear transformation, where \code{"mean"} relies on the intercept and \code{"linear"} (default) the intercept and slope
}
  \item{method}{
string specifying the equating method, where \code{"none"} (default) indicates observed score equating under the random groups design, and both \code{"tucker"} and \code{"levine"} indicate the respective methods under the nonequivalent groups design
}
  \item{xv, yv}{
vectors of anchor test total scores for form X and Y examinees
}
  \item{vscale}{
the score scale for the anchor test, as a vector
}
  \item{w}{
the weight applied to form X when calculating the synthetic population
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
As noted above, for the random groups design \code{x} and \code{y} can be supplied as either total scores, where the length of each will equal the number of examinees, or as counts, where the length of each will equal the length of the score scale, with each value indicating the total number of examinees obtaining the corresponding score. With the nonequivalent groups design, only total scores may be used.

Two equating types (mean and linear) and two equating methods (Tucker and Levine) are currently supported. Varying the \code{type} and \code{method} arguments results in four possible equating scenarios. Whenever \code{"tucker"} or \code{"levine"} is specified, nonequivalent groups with anchor test design is assumed and the anchor test arguments \code{xv}, \code{yv}, \code{vscale} and \code{w} are also required (see \code{\link{synthetic}} for details).
}
\value{
Returns a list including the following components when \code{verbose = TRUE} (only the first is returned when \code{verbose = FALSE}):
  \item{yx}{form Y equivalents of the form X scores. For the random groups design, when \code{verbose = TRUE} a column of standard errors is included}
  \item{coefficients}{the slope and intercept used in the linear conversion}
  \item{synthstats}{matrix of means and standard deviations for the anchor and synthetic distributions}
  \item{anchortab}{frequency table for the anchor test}
}
\references{
Kolen, M. J., & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\note{}
\seealso{\code{\link{equate}}, \code{\link{synthetic}}}
\examples{
# see the function 'equate' for more examples
x <- ACTmath[,2]
y <- ACTmath[,3]
scale <- ACTmath[,1]
cbind(scale,equate.ln(x,y,scale,type="m")$yx,
  equate.ln(x,y,scale,type="l")$yx)
}
\keyword{methods}