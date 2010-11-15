\name{equate}
\alias{equate}
\alias{equate.ln}
\alias{equate.eq}
\alias{equate.ca}
\alias{equipercentile}
\alias{synthetic}
\alias{se.eq}
\alias{se.ln}
\alias{print.equate}
\title{Statistical Equating}
\description{
This function equates test scores under the random groups and nonequivalent groups with anchor test designs.
}
\usage{
equate(x, y, type, method = NA, name = NULL, bootse = FALSE, ...)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y. These may be univariate or bivariate distributions (see below for details). Alternatively, if \code{y} is an object of class \dQuote{\code{\link{equate}}} (the output of a previous equating) and \code{x} a vector of scores, the equating function from \code{y} is used to convert scores in \code{x} to the Y scale.
}
  \item{type}{
the type of equating, as the string \code{"mean"}, \code{"linear"}, \code{"equipercentile"}, \code{"circle-arc"}, or abbreviations of these
}
  \item{method}{
string specifying the equating method, where \code{NA} (default) indicates equating under the random groups design, and \code{"nominal weights"}, \code{"tucker"}, \code{"levine"}, \code{"frequency estimation"}, \code{"braun/holland"}, and \code{"chained"} (or abbreviations of these) indicate the respective methods under the nonequivalent groups design
}
  \item{name}{
An optional name, used to label the output
}
  \item{bootse}{
logical indicating whether or not bootstrap standard errors of equating should be returned (default is \code{FALSE})
}
  \item{\dots}{
further arguments passed to or from other methods (see below for details)
}
}
\details{
For the random groups design \code{x} and \code{y} are separate frequency tables where each row contains a score scale point and the total number of examinees obtaining that score. For the nonequivalent groups design each row contains a score scale point on \code{x} or \code{y}, a score scale point on the anchor test, and the total number of examinees obtaining that score combination (see \code{\link{freqtab}} for details). When equating output is supplied for \code{y}, no other arguments are required and scores from \code{x} are converted directly to the scale indicated in \code{y}. In this case \code{x} must be a vector of scores.

Four equating types are currently supported: mean, linear, and equipercentile, as described by Kolen and Brennan (2004), and circle-arc equating, as described by Livingston and Kim (2009). The equating design is implied by the \code{method} argument, where \code{NA} (default) indicates that no method is needed (because examinees taking forms X and Y are assumed to be equivalent). The nominal weights, Tucker, Levine observed score, Levine true score, frequency estimation, Braun/Holland, and chained equating methods are supported for the nonequivalent groups with anchor test design. All but the Levine true score and chained method rely on a \dQuote{synthetic} distribution of scores (Braun and Holland, 1982), a weighted combination of \code{x} and \code{y}.

Depending on the equating method, the following additional arguments may be required:
\describe{
  \item{\code{internal}}{
logical indicating whether or not the anchor item scores are included in the total scores. This applies only to the Levine method, as all other methods assume an internal anchor test. Default is \code{TRUE}
}
  \item{\code{Ky}}{
integer, typically indicating the number of items on form Y or the maximum possible score on Y (default). This value is used in equipercentile equating for the extremes of the scale (see below)
}
  \item{\code{smoothmeth}}{
string indicating one of four smoothing methods to be used in equipercentile equating: \code{"none"} (default), \code{"bump"}, \code{"average"}, and \code{"loglin"} (see below)
}
  \item{\code{w}}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w}) when estimating the synthetic population. When set to -1, proportional weights are calculated for X and Y based on sample size
}
  \item{\code{lts}}{
logical indicating whether or not to use levine true score (hence \dQuote{lts}) equating. Default is \code{FALSE}
}
  \item{\code{lowp}, \code{highp}}{
two vectors, each of length 2, specifying the coordinates for the low and high points used in circle-arc equating. \code{lowp} defaults to the minimums of the X and Y scales. \code{highp} defaults to the number of items, i.e., the maximum possible scores
}
  \item{\code{midp}}{
string used to specify the type of chained linear equating used to obtain the midpoint in circle-arc equating, whether \code{"mean"} (default) or \code{"linear"}
}
}
For equipercentile equating under the random groups design, three smoothing options are available: \code{"bump"} and \code{"average"} (see \code{\link{freqbump}} and \code{\link{freqavg}}) require the additional argument \code{jmin}, and loglinear smoothing (\code{"loglin"}) requires the score functions \code{xscorefun} and \code{yscorefun} (see \code{\link{loglinear}}). For frequency estimation and chained methods, only smoothing methods \code{"freqbump"} and \code{"loglin"} are supported.

The parameter \code{Ky} is used to obtain form Y equivalents of form X scores with percentile ranks of 100. Typically \code{Ky} is set to be the number of score points in the form Y scale, which assumes that scores are integers ranging from 1 (or 0) to the total number of items, and that each item is scored correct/incorrect. Scores on other scales (such as scales which include negative values, or which exclude zero) may also be used. In such cases \code{Ky} can be set to the highest possible score on form Y, or alternatively the highest observed score on Y. Specifying \code{Ky = "maxobs"} sets Ky to the maximum observed score on Y. This may be useful for calculating bootstrap standard errors, where Ky would depend on the particular bootstrap sample.

Standard errors are returned for the random groups design only, using equations derived by Lord (1982, p. 168) for equipercentile equating and Braun and Holland (1982) for linear. Bootstrap standard errors are estimated using a default of \code{reps=100} replications, sampling the maximum amount from each score distribution (controlled by the arguments \code{xn} and \code{yn}; see \code{\link{se.boot}}).
}
\value{
An object of class \dQuote{\code{\link{equate}}} is returned, listing the following components, some of which are conditional on the equating type, method, and smoothing:
  \item{type}{equating type}
  \item{method}{equating method}
  \item{design}{equating design, inferred from the method}
  \item{stats}{descriptive statistics for X, Y, and YX (equated scores) distributions}
  \item{freqtab}{table containing the score scale, frequencies, and the cumulative distribution functions for forms X and Y}
  \item{yxtab}{frequency table, of class \dQuote{freqtab}, for equated scores}
  \item{concordance}{conversion table containing scores on X with their form Y equivalents. For the random groups design a column of standard errors is included}
  \item{coefficients}{conversion coefficients intercept and slope. Center points of circle and radius are also included for circle-arc equating}
  \item{points}{X and Y coordinates for low, middle, and high points used in circle-arc equating}
  \item{synthstats}{means and standard deviations for the synthetic distributions}
  \item{w}{weight applied to X in synthetic estimation}
  \item{synthtab}{frequency table for the synthetic distribution, used in frequency estimation equipercentile equating}
  \item{anchorstats}{descriptive statistics for the anchor distributions}
  \item{anchortab}{frequency table for the anchor distributions}
  \item{smoothmethod}{smoothing method}
  \item{smoothout}{a list of smoothing output, the length of which depends on the method chosen (for loglinear smoothing this list can be extremely long)}
  \item{bootsee}{vector of bootstrap standard errors of equating}
}
\references{
Braun, H. I., and Holland, P. W. (1982). Observed-score test equating: A mathematical analysis of some ETS equating procedures. In P. W. Holland and D. B. Rubin (Eds.), \emph{Test Equating} (pp. 9--49). New York: Academic.

Kolen, M. J., and Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.

Livingston, S. A., and Kim, S. (2009). The circle-arc method for equating in small samples, \emph{Journal of Educational Measurement, 46}, 330--343.

Lord, F. (1982). The standard error of equipercentile equating. \emph{Journal of Educational Statistics, 7}, 165--174.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqbump}}, \code{\link{freqavg}}, \code{\link{loglinear}}, \code{\link{se.boot}}
}
\examples{
# see vignette("equatevignette") for description of methods
# and additional examples

# random groups equating for (1) mean, (2) linear, 
# (3) equipercentile, and (4) equipercentile with 
# loglinear smoothing:
rscale <- ACTmath[, 1]
rx <- as.freqtab(rscale, ACTmath[, 2])
ry <- as.freqtab(rscale, ACTmath[, 3])
set.seed(2007)

req1 <- equate(rx, ry, type = "m", bootse = TRUE, reps = 50)
req2 <- equate(rx, ry, type = "l", bootse = TRUE, reps = 50)
req3 <- equate(rx, ry, type = "e", bootse = TRUE, reps = 50)
req4 <- equate(rx, ry, type = "e", bootse = TRUE, reps = 50,
  smooth = "loglin", degree = 3)

# compare equated scores:
cbind(rscale, mean = req1$conc[, 2], linear = req2$conc[, 2],
  equip = req3$conc[, 2], equipS = req4$conc[, 2])

# compare boostrap standard errors:
cbind(rscale, linear = req2$bootsee, equip = req3$bootsee,
  equipS = req4$bootsee)


# nonequivalent groups design for (1) Tucker linear,
# (2) frequency estimation , and (3) Braun/Holland linear
nscale <- 0:36
vscale <- 0:12
nx <- freqtab(nscale, KBneat$x[, 1], vscale, KBneat$x[, 2])
ny <- freqtab(nscale, KBneat$y[, 1], vscale, KBneat$y[, 2])

neq1 <- equate(nx, ny, type = "Linear", method = "Tuck", w = 1)
neq2 <- equate(nx, ny, type = "equip", method = "freq", w = 1)
neq3 <- equate(nx, ny, type = "linear", method = "braun", w = 1)

# compare equated scores:
cbind(nscale, Tucker = neq1$conc[, 2], FE = neq2$conc[, 2],
  Braun = neq3$conc[, 2])
}
\keyword{methods}
