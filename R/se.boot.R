se.boot <- function(x, y, xn = sum(x[, ncol(x)]),
  yn = sum(y[, ncol(y)]), reps = 100, eqfun,
  returnboots = FALSE, ...) {

  nc <- ncol(x)
  xscale <- unique(x[, 1])
  xprob <- x[, nc]/sum(x[, nc])
  yprob <- y[, nc]/sum(y[, nc])
  eqfun <- match.fun(eqfun)
  eqmat <- matrix(nrow = length(xscale), ncol = reps)
  if(nc == 3) {
    index <- 1:nrow(x)
    vscale <- unique(x[, 2])
    for(i in 1:reps) {
      xi <- sample(index, xn, replace = TRUE, prob = xprob)
      xtemp <- freqtab(xscale, x[xi, 1], vscale, x[xi, 2])
      yi <- sample(index, yn, replace = TRUE, prob = yprob)
      ytemp <- freqtab(xscale, y[yi, 1], vscale, y[yi, 2])
      eqmat[, i] <- eqfun(xtemp, ytemp, ...)
    }
  }
  else {
    for(i in 1:reps) {
      xtemp <- freqtab(xscale, sample(x[, 1], xn,
        replace = TRUE, prob = xprob))
      ytemp <- freqtab(xscale, sample(y[, 1], yn,
        replace = TRUE, prob = yprob))
      eqmat[, i] <- eqfun(xtemp, ytemp, ...)
    }
  }
  if(returnboots)
    return(eqmat)
  else return(apply(eqmat, 1, sd))
}
