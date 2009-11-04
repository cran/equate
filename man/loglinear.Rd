\name{loglinear}
\Rdversion{1.1}
\alias{loglinear}
\title{Loglinear Smoothing}
\description{
This function smooths univariate and bivariate score distributions via polynomial loglinear modeling (Holland & Thayer, 2000; Moses & von Davier, 2006).
}
\usage{
loglinear(x, scorefun, degree, raw = TRUE, convergecrit = .0001, ...)
}
\arguments{
  \item{x}{
score distribution of class \dQuote{\code{\link{freqtab}}} for form X. Under the random groups design (i.e., no anchor test) \code{x} will contain the score scale in column 1, and the number of examinees obtaining each score in column 2. For the nonequivalent groups design a bivariate frequency table is used, where columns 1 and 2 include all score combinations for the total and anchor test score scales, and column 3 contains the number of examinees obtaining each combination (see \code{\link{freqtab}} for details)
}
  \item{scorefun}{
matrix of score functions, where each column represents a transformation of the score scale (or the crossed score scales, in the bivariate case)
}
  \item{degree}{
integer indicating a maximum polynomial transformation to be computed (passed to \code{\link{poly}}; ignored if \code{scorefun} is provided)
}
  \item{raw}{
logical. If \code{TRUE} (default), raw polynomials will be used, if \code{FALSE}, orthogonal (passed to \code{\link{poly}})
}
  \item{convergecrit}{
convergence criteria used in maximum likelihood estimation (default is .0001)
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
Loglinear smoothing is a flexible procedure for reducing irregularities in a raw score distribution. \code{loglinear} fits a polynomial loglinear model to a distribution of scores, where the degree of each polynomial term determines the specific moment of the raw distribution that is preserved in the fitted distribution (see below for examples). \code{scorefun} must contain at least one score function of the scale score values. While there is no explicit limit on the number of columns in \code{scorefun}, models with more than ten may not converge on a solution (the complexity of the model must also be taken into account).

In the univariate case, specifying \code{degree} is an alternative to \code{scorefun}. It takes a maximum polynomial degree and constructs the score functions accordingly. For example, \code{degree=3} will result in a model with three terms: the score scale raised to the first, second, and third powers.

Maximum likelihood estimates are obtained using a Newton-Raphson algorithm, with slightly smoothed frequencies (all nonzero) as the basis for starting values. Calculating standard errors for these estimates requires matrix inversion, which for complex models may not be possible. In this case the standard errors will be omitted. The tolerance level for detecting singularity may be modified with the argument \code{tol}, which is passed to \code{solve()}.

For a detailed description of the estimation procedures, including examples, see Holland and Thayer, 1987 and 2000. For a more recent discussion, including the macro after which the \code{loglinear} function is modeled, see Moses and von Davier, 2006.
}
\value{
Returns a list including the following components:
  \item{modelfit}{table of model fit statistics: likelihood ratio chi-square, Pearson chi-square, Freeman-Tukey chi-square, AIC, and CAIC}
  \item{rawbetas}{two-column matrix of raw maximum likelihood estimates for the beta coefficients and corresponding standard errors}
  \item{alpha}{normalizing constant}
  \item{iterations}{number of iterations reached before convergence}
  \item{fitted.values}{vector of estimated frequencies}
  \item{residuals}{vector of residuals}
  \item{cmatrix}{the \dQuote{C matrix}, a factorization of the covariance matrix of the fitted values}
  \item{scorefun}{matrix of score functions}
}
\references{
Holland, P. W., & Thayer, D. T. (1987). \emph{Notes on the use of log-linear models for fitting discrete probability distributions} (PSR Technical Rep. No. 87-79; ETS RR-87-31). Princeton, NJ: ETS.

Holland, P. W., & Thayer, D. T. (2000). Univariate and bivariate loglinear models for discrete test score distributions. \emph{Journal of Educational and Behavioral Statistics, 25}, 133-183.

Moses, T., & von Davier, A. A. (2008). \emph{A SAS macro for loglinear smoothing: Applications and implications} (ETS Research Rep. No. RR-08-59). Princeton, NJ: ETS.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{\code{\link{glm}}, \code{\link{loglin}}
}
\examples{
set.seed(2010)
x <- round(rnorm(1000,100,15))
xscale <- 50:150

# smooth x preserving first 3 moments:
xlog1 <- loglinear(freqtab(x,xscale),degree=3)
xtab <- freqtab(x,xscale)
cbind(xtab,xlog1$fit)

par(mfrow=c(2,1))
plot(xscale,xtab[,2],type="h",ylab="count",
  main="X raw")
plot(xscale,xlog1$fit,type="h",ylab="count",
  main="X smooth")

# add "teeth" and "gaps" to x:
teeth <- c(.5,rep(c(1,1,1,1,.5),20))
xt <- xtab[,2]*teeth
cbind(xtab,xt)
xttab <- freqtab(xt,xscale,add=TRUE)
xlog2 <- loglinear(xttab,degree=3)
cbind(xscale,xt,xlog2$fit)

# smooth xt using score functions that preserve 
# the teeth structure (also 3 moments):
teeth2 <- c(1,rep(c(0,0,0,0,1),20))
xt.fun <- cbind(xscale,xscale^2,xscale^3)
xt.fun <- cbind(xt.fun,xt.fun*teeth2)
xlog3 <- loglinear(xttab,xt.fun)
cbind(xscale,xt,xlog3$fit)

par(mfrow=c(2,1))
plot(xscale,xt,type="h",ylab="count",
  main="X teeth raw")
plot(xscale,xlog3$fit,type="h",ylab="count",
  main="X teeth smooth")

# bivariate example, preserving first 3 moments of y
# and v (anchor) each, the covariance of y and v, and
# 3 additional degrees of dependence in y and v:
yv <- KBneat$y
yscale <- 0:36
vscale <- 0:12
yvtab <- freqtab(yv[,1],yscale,yv[,2],vscale)
Y <- yvtab[,1]
V <- yvtab[,2]
scorefun <- cbind(Y,Y^2,Y^3,V,V^2,V^3,V*Y,V*Y*Y,V*V*Y,V*V*Y*Y)
loglinear(yvtab,scorefun)$model

# replicate Moses and von Davier, 2006, univariate example:
uv <- c(0,4,11,16,18,34,63,89,87,129,124,154,125,
  131,109,98,89,66,54,37,17)
loglinear(freqtab(uv,0:20,add=TRUE),degree=3)
}
\keyword{smooth}