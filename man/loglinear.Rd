\name{loglinear}
\alias{loglinear}
\title{
Polynomial Loglinear Smoothing
}
\description{
This function calls on \code{\link{glm}} to fit one or more loglinear smoothing models to a frequency distribution. Model fit can also be compared for nested models.
}
\usage{
loglinear(x, scorefun, degree, raw = TRUE, verbose = FALSE,
	compare = FALSE, stepup = compare, showWarnings = FALSE, ...)
}
\arguments{
  \item{x}{
score distribution of class \dQuote{\code{\link{freqtab}}}. With an equivalent groups design (i.e., no anchor test), \code{x} will contain the score scale in column 1 and the number of examinees obtaining each score in column 2. With a nonequivalent groups design, a bivariate frequency table is used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details).
}
  \item{scorefun}{
matrix of score functions, where each column includes a transformation of the score scale (or the crossed score scales, in the bivariate case). If missing, \code{degree} will be used to construct polynomial score functions.
}
  \item{degree}{
integer indicating a maximum polynomial score transformation to be computed (passed to \code{\link{poly}}; ignored if \code{scorefun} is provided).
}
  \item{raw}{
logical, passed to \code{\link{poly}}, indicating whether raw polynomials (\code{TRUE}, default) or orthogonal polynomials (\code{FALSE}) will be used.
}
  \item{verbose}{
logical, with default \code{TRUE}, indicating whether or not \code{glm} output in addition to the fitted values should be returned.
}
  \item{compare}{
logical, with default \code{FALSE}, indicating whether or not fit for nested models should be compared. If \code{TRUE}, \code{stepup} is also set to \code{TRUE} and only results from the model fit comparison are returned.
}
  \item{stepup}{
logical, with default \code{FALSE}, indicating whether or not all nested models should also be run.
}
  \item{showWarnings}{
logical, with default \code{FALSE}, indicating whether or not warnings from \code{glm} should be shown.
}
  \item{\dots}{
further arguments passed to \code{glm}.
}
}
\details{
Loglinear smoothing is a flexible procedure for reducing irregularities in a frequency distribution prior to equating, where the degree of each polynomial term determines the specific moment of the raw distribution that is preserved in the fitted distribution (see below for examples). The \code{loglinear} function is a simple wrapper for \code{\link{glm}}, and is used to simplify the creation of polynomial score functions and the fitting and comparing of multiple loglinear models. \code{scorefun} must contain at least one score function of the scale score values.

Specifying \code{degree} is an alternative to supplying \code{scorefun}. \code{degree} takes a maximum polynomial degree and constructs the score functions accordingly. For example, \code{degree = 3} will result in a model with three terms: the score scale raised to the first, second, and third powers, preserving the mean, variance, and skew of the original distribution. In the bivariate case, powers 1 through 3 would be included for each variable.

\code{stepup} is used to run nested models based on subsets of the columns in \code{scorefun}. Output will correspond to models based on columns 1 and 2, 1 through 3, 1 through 4, etc.

\code{compare} returns output from \code{anova}, comparing model fit for all the models run with \code{stepup = TRUE}.

For additional examples, see Holland and Thayer, 1987 and 2000.
}
\value{
Returns either an anova fit table (when \code{compare = TRUE}), a vector or matrix of fitted values (when \code{verbose = FALSE}), or complete output from \code{glm} (when \code{verbose = TRUE}) including for nested models (when \code{stepup = TRUE}).
}
\references{
Holland, P. W., and Thayer, D. T. (1987). \emph{Notes on the use of log-linear models for fitting discrete probability distributions} (PSR Technical Rep. No. 87-79; ETS RR-87-31). Princeton, NJ: ETS.

Holland, P. W., and Thayer, D. T. (2000). Univariate and bivariate loglinear models for discrete test score distributions. \emph{Journal of Educational and Behavioral Statistics, 25}, 133--183.

Moses, T. P., and von Davier, A. A. (2006). \emph{A SAS macro for loglinear smoothing: Applications and implications} (ETS Research Rep. No. RR-06-05). Princeton, NJ: ETS.

Wang, T. (2009). Standard errors of equating for the percentile rank-based equipercentile equating with log-linear presmoothing. \emph{Journal of Educational and Behavioral Statistics, 34}, 7--23.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{glm}}, \code{anova.glm}, \code{\link{loglin}}
}
\examples{
set.seed(2010)
x <- round(rnorm(1000, 100, 15))
xscale <- 50:150

# smooth x preserving first 3 moments:
xtab <- freqtab(x, xscale = xscale)
xlog1 <- loglinear(xtab, degree = 3)
plot(xtab, y = xlog1)
lines(xtab[, 1], xlog1)

# add "teeth" and "gaps" to x:
teeth <- c(.5, rep(c(1, 1, 1, 1, .5), 20))
xt <- xtab[, 2] * teeth
cbind(xtab, xt)
xttab <- as.freqtab(cbind(xscale, xt))
xlog2 <- loglinear(xttab, degree = 3)
cbind(xscale, xt, xlog2)

# smooth xt using score functions that preserve 
# the teeth structure (also 3 moments):
teeth2 <- c(1, rep(c(0, 0, 0, 0, 1), 20))
xt.fun <- cbind(xscale, xscale^2, xscale^3)
xt.fun <- cbind(xt.fun, xt.fun * teeth2)
xlog3 <- loglinear(xttab, xt.fun)
cbind(xscale, xt, xlog3)

par(mfrow = c(2, 1))
plot(xscale, xt, type = "h", ylab = "count",
	main = "X teeth raw")
plot(xscale, xlog3, type = "h", ylab = "count",
	main = "X teeth smooth")

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
scorefun <- cbind(Y, Y^2, Y^3, V, V^2, V^3, V * Y)
wang09 <- equate(xvtab, yvtab, type = "equip",
	method = "chained", smooth = "loglin", xscorefun = scorefun, 
	yscorefun = scorefun)
wang09$concordance

# replicate Moses and von Davier, 2006, univariate example:
uv <- c(0, 4, 11, 16, 18, 34, 63, 89, 87, 129, 124,
	154, 125, 131, 109, 98, 89, 66, 54, 37, 17)
loglinear(as.freqtab(cbind(0:20, uv)), degree = 3)
}
\keyword{smooth}
\keyword{models}
