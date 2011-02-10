equate.eq <- function(x, y, method = NA, ident, Ky = max(y[, 1]),
  w = -1, smoothmeth = "none", jmin, xscorefun, yscorefun,
  verbose = FALSE, ...) {

  nc <- ncol(x)
  xscale <- unique(x[, 1])
  yscale <- unique(y[, 1])
  xcount <- as.vector(tapply(x[, nc], x[, 1], sum))
  ycount <- as.vector(tapply(y[, nc], y[, 1], sum))
  xtab <- x
  ytab <- y

  smoothmeth <- match.arg(tolower(smoothmeth), c("none",
    "bump", "average", "loglin"))
  if(smoothmeth == "bump") {
    if(missing(jmin))
      jmin <- 10^-6
    xsmooth <- freqbump(x, jmin) * sum(x[, nc])
    ysmooth <- freqbump(y, jmin) * sum(y[, nc])
    xtab[, nc] <- xsmooth
    ytab[, nc] <- ysmooth
  }
  else if(smoothmeth == "average" & nc == 2) {
    if(missing(jmin))
      jmin <- 1
    xsmooth <- freqavg(x, jmin)
    ysmooth <- freqavg(y, jmin)
    xtab[, nc] <- xsmooth
    ytab[, nc] <- ysmooth
  }
  else if(smoothmeth == "loglin") {
    xsmooth <- loglinear(x, xscorefun, ...)$fitted
    ysmooth <- loglinear(y, yscorefun, ...)$fitted
    xtab[, nc] <- xsmooth
    ytab[, nc] <- ysmooth
  }

  method <- match.arg(tolower(method),
    c(NA, "frequency estimation", "chained"))
  if(is.na(method))
    yx <- equipercentile(xtab, ytab, Ky)
  else if(method == "frequency estimation") {
    stabs <- synthetic(xtab, ytab, w, method)
    xtab <- as.freqtab(xscale, stabs$synthtab[, 2])
    ytab <- as.freqtab(yscale, stabs$synthtab[, 3])
    yx <- equipercentile(xtab, ytab, Ky)
  }
  else if(method == "chained") {
    xvtab <- as.freqtab(unique(xtab[, 2]),
      tapply(xtab[, 3], xtab[, 2], sum))
    xtab <- as.freqtab(xscale, tapply(xtab[, 3], xtab[, 1], sum))
    yvtab <- as.freqtab(unique(ytab[, 2]),
      tapply(ytab[, 3], ytab[, 2], sum))
    ytab <- as.freqtab(yscale, tapply(ytab[, 3], ytab[, 1], sum))
    vx <- equipercentile(xtab, xvtab)$yx
    pvyx <- px(vx, yvtab)
    yx <- equipercentile(pvyx, ytab)
  }

  out <- list(yx = (1 - ident) * yx$yx + ident * xscale)

  if(verbose) {
    if(is.na(method))
      out$se <- yx$se
    out$stats <- rbind(x = c(descript(x)), y = c(descript(y)),
      yx = c(descript(as.freqtab(out$yx, xcount))))
    colnames(out$stats) <- c("mean", "sd", "skew", "kurt", "n")
    out$freqtab <- cbind(scale = xscale, fx = fx(xcount),
      fy = fx(ycount), xcount = xcount, ycount = ycount)
    if(!is.na(method)) {
      out$w <- w
      out$anchorstats <- rbind(descript(x[, -1]),
        descript(y[, -1]))
      rownames(out$anchorstats) <- c("xv", "yv")
      out$anchortab <- cbind(scale = unique(x[, 2]),
        xvcount = tapply(x[, 3], x[, 2], sum),
        yvcount = tapply(y[, 3], y[, 2],sum))
      if(method == "frequency estimation")
        out <- c(out, stabs)
      if(method == "chained")
        out$chaintab <- cbind(vxx = vx, pvyx = pvyx)
    }
    if(smoothmeth != "none") {
      out$smoothmethod <- smoothmeth
      out$smoothout <- list(x = xsmooth, y = ysmooth)
      if(!is.na(method)) {
        out$anchortab <- cbind(out$anchortab,
          xvfitted = tapply(xsmooth, x[, 2], sum),
          yvfitted = tapply(ysmooth, y[, 2], sum))
        xsmooth <- tapply(xsmooth, x[, 1], sum)
        ysmooth <- tapply(ysmooth, y[, 1], sum)
      }
      out$freqtab <- cbind(out$freqtab, xfitted = xsmooth,
        yfitted = ysmooth)
    }
  }
  else out <- out$yx
  return(out)
}
