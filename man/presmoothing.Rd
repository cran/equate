\name{Presmoothing}
\alias{presmoothing}
\alias{loglinear}
\alias{freqbump}
\alias{freqavg}
\title{
Frequency Distribution Presmoothing
}
\description{
These functions are used to smooth frequency distributions.
}
\usage{
presmoothing(x, smoothmethod = c("none", "average", "bump",
	"loglinear"), jmin,	asfreqtab = TRUE, ...)

loglinear(x, scorefun, degree = 4, xdegree = 1, raw = TRUE,
	asfreqtab = TRUE, verbose = FALSE, compare = FALSE,
	stepup = compare, showWarnings = FALSE, ...)
}
\arguments{
  \item{x}{
score distribution of class \dQuote{\code{\link{freqtab}}}. In the univariate case, \code{x} will contain the score scale in column 1 and the number of examinees obtaining each score in column 2. With a single group or nonequivalent groups design, a bivariate frequency table is used, where columns 1 and 2 include all score combinations for the two score scales, and column 3 contains the number of examinees obtaining each combination.
}
  \item{smoothmethod}{
character string indicating the smoothing method to be used by \code{presmoothing}. \code{"none"} returns unsmoothed frequencies, \code{"bump"} adds a small frequency to each score value, \code{"average"} imputes small frequencies with average values, and \code{"loglinear"} fits loglinear models. See below for details.
}
  \item{jmin}{
for \code{smoothmethod = "average"}, the minimum frequency, as an integer, below which frequencies will be replaced (default is 1).
for \code{smoothmethod = "bump"}, the value to be added to each score point (as a probability, with default 1e-6).
}
  \item{asfreqtab}{
logical, with default \code{TRUE}, indicating whether or not a frequency table should be returned. For \code{smoothmethod = "average"} and \code{smoothmethod = "bump"}, the alternative is a vector of frequencies. For \code{loglinear}, there are other alternatives.
}
  \item{scorefun}{
matrix of score functions, where each column includes a transformation of the score scale (or the crossed score scales, in the bivariate case). If missing, \code{degree} and \code{xdegree} will be used to construct polynomial score functions.
}
  \item{degree, xdegree}{
integers indicating the maximum polynomial score transformation to be computed. Defaults are \code{4} and \code{1}. These are ignored if \code{scorefun} is provided; otherwise, they are passed to \code{\link{poly}}.
}
  \item{raw}{
logical, passed to \code{\link{poly}}, indicating whether raw polynomials (\code{TRUE}, default) or orthogonal polynomials (\code{FALSE}) will be used.
}
  \item{verbose}{
logical, with default \code{FALSE}, indicating whether or not full \code{glm} output should be returned.
}
  \item{compare}{
logical, with default \code{FALSE}, indicating whether or not fit for nested models should be compared. If \code{TRUE}, \code{stepup} is also set to \code{TRUE} and only results from the model fit comparison are returned, that is, \code{verbose} is ignored.
}
  \item{stepup}{
logical indicating whether or not all nested models should also be run.
}
  \item{showWarnings}{
logical, with default \code{FALSE}, indicating whether or not warnings from \code{glm} should be shown.
}
  \item{\dots}{
further arguments passed to other methods. For \code{presmoothing}, these are passed to \code{loglinear} and include those listed above.
}
}
\details{
Loglinear smoothing is a flexible procedure for reducing irregularities in a frequency distribution prior to equating, where the degree of each polynomial term determines the specific moment of the observed distribution that is preserved in the fitted distribution (see below for examples). The \code{loglinear} function is a wrapper for \code{\link{glm}}, and is used to simplify the creation of polynomial score functions and the fitting and comparing of multiple loglinear models.

\code{scorefun}, if supplied, must contain at least one score function of the scale score values. Specifying \code{degree} and \code{xdegree} is an alternative to supplying \code{scorefun}. These arguments take the maximum polynomial degree for the univariate (\code{degree}) and bivariate (\code{xdegree}) distributions, and constructs the score functions accordingly. For example, \code{degree = 3} will result in a model with three terms: the score scale raised to the first, second, and third powers, preserving the mean, variance, and skew of the original distribution. In the bivariate case, powers 1 through 3 would be included for both variables. And \code{xdegree = 1} would add a polynomial for the interaction of the two score scales, preserving their covariance.

\code{stepup} is used to run nested models based on subsets of the columns in \code{scorefun}. Output will correspond to models based on columns 1 and 2, 1 through 3, 1 through 4, to 1 through \code{ncol(scorefun)}.

\code{compare} returns output as an \code{anova} table, comparing model fit for all the models run with \code{stepup = TRUE}.

The remaining smoothing methods make adjustments to scores with low or zero frequencies. \code{smoothmethod = "bump"} adds the proportion \code{jmin} to each score point and then adjusts the probabilities to sum to 1. \code{smoothmethod = "average"} replaces frequencies falling below the minimum \code{jmin} with averages of adjacent values.
}
\value{
When \code{smoothmethod = "average"} or \code{smoothmethod = "bump"}, either a smoothed frequency vector or table is returned. Otherwise, \code{loglinear} returns the following:
\item{}{when \code{compare = TRUE}, an anova table for model fit}
\item{}{when \code{asfreqtab = TRUE}, a smoothed frequency table}
\item{}{when \code{verbose = TRUE}, full \code{glm} output, for all nested models when \code{stepup = TRUE}}
\item{}{when \code{stepup = TRUE} and \code{verbose = FALSE}, a \code{data.frame} of fitted frequencies, with one column per model}
}
\references{
Holland, P. W., and Thayer, D. T. (1987). \emph{Notes on the use of log-linear models for fitting discrete probability distributions} (PSR Technical Rep. No. 87-79; ETS RR-87-31). Princeton, NJ: ETS.

Holland, P. W., and Thayer, D. T. (2000). Univariate and bivariate loglinear models for discrete test score distributions. \emph{Journal of Educational and Behavioral Statistics, 25}, 133--183.

Moses, T. P., and von Davier, A. A. (2006). \emph{A SAS macro for loglinear smoothing: Applications and implications} (ETS Research Rep. No. RR-06-05). Princeton, NJ: ETS.

Moses, T., and Holland, P. W. (2008). \emph{Notes on a general framework for observed score equating} (ETS Research Rep. No. RR-08-59). Princeton, NJ: ETS.

Wang, T. (2009). Standard errors of equating for the percentile rank-based equipercentile equating with log-linear presmoothing. \emph{Journal of Educational and Behavioral Statistics, 34}, 7--23.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{glm}}, \code{\link{loglin}}
}
\examples{
set.seed(2010)
x <- round(rnorm(1000, 100, 15))
xscale <- 50:150
xtab <- freqtab(x, xscale = xscale)

# adjust frequencies
plot(xtab, y = cbind(average = freqavg(xtab),
  bump = freqbump(xtab)))

# smooth x preserving first 3 moments:
xlog1 <- loglinear(xtab, degree = 3,
  asfreqtab = FALSE)
plot(xtab, y = xlog1, legendtext = "degree = 3")

# add "teeth" and "gaps" to x:
teeth <- c(.5, rep(c(1, 1, 1, 1, .5), 20))
xt <- xtab[, 2]*teeth
xttab <- as.freqtab(cbind(xscale, xt))
xlog2 <- loglinear(xttab, degree = 3)

# smooth xt using score functions that preserve 
# the teeth structure (also 3 moments):
teeth2 <- c(1, rep(c(0, 0, 0, 0, 1), 20))
xt.fun <- cbind(xscale, xscale^2, xscale^3)
xt.fun <- cbind(xt.fun, teeth2, xt.fun * teeth2)
xlog3 <- loglinear(xttab, xt.fun)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
plot(xttab, main = "unsmoothed", ylim = c(0, 30))
plot(xlog2, main = "ignoring teeth", ylim = c(0, 30))
plot(xlog3, main = "preserving teeth", ylim = c(0, 30))
par(op)

# bivariate example, preserving first 3 moments of total
# and v (anchor) each of x and y, and the covariance
# between anchor and total
# see equated scores in Wang (2009), Table 4
xvtab <- freqtab(KBneat$x[, 1], KBneat$x[, 2],
  xscale = 0:36, vscale = 0:12)
yvtab <- freqtab(KBneat$y[, 1], KBneat$y[, 2],
  xscale = 0:36, vscale = 0:12)
Y <- yvtab[, 1]
V <- yvtab[, 2]
scorefun <- cbind(Y, Y^2, Y^3, V, V^2, V^3, V*Y)
wang09 <- equate(xvtab, yvtab, type = "equip",
  method = "chained", smooth = "loglin",
  scorefun = scorefun)
wang09$concordance

# replicate Moses and von Davier, 2006, univariate example:
uv <- c(0, 4, 11, 16, 18, 34, 63, 89, 87, 129, 124,
  154, 125, 131, 109, 98, 89, 66, 54, 37, 17)
loglinear(as.freqtab(cbind(0:20, uv)), degree = 3)
}
\keyword{smooth}
\keyword{models}
