equate.eq <- function(x, y, method = NA, Ky = max(y[, 1]),
  w = -1, smoothmeth = "none", jmin, xscorefun, yscorefun,
  verbose = FALSE, ...) {

  nc <- ncol(x)
  xscale <- unique(x[, 1])
  yscale <- unique(y[, 1])
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
    xtab[, 2] <- xsmooth
    ytab[, 2] <- ysmooth
  }
  else if(smoothmeth == "loglin") {
    xsmooth <- loglinear(x, xscorefun, ...)
    ysmooth <- loglinear(y, yscorefun, ...)
    xtab[, nc] <- xsmooth$fitted
    ytab[, nc] <- ysmooth$fitted
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
    vx <- equipercentile(xtab, xvtab)[, 1]
    pvyx <- px(vx, yvtab)
    yx <- cbind(equipercentile(pvyx, ytab))
  }

  if(verbose) {
    if(is.na(method))
      out <- list(yx=yx)
    else {
      out <- list(yx = yx[, 1])
      if(method == "frequency estimation")
        out <- c(out, stabs)
      if(method == "chained")
        out$chaintab <- cbind(vxx = vx, pvyx = pvyx)
      out$anchorstats <- rbind(descript(x[, -1]),
        descript(y[, -1]))
      rownames(out$anchorstats) <- c("xv", "yv")
      out$anchortab <- cbind(scale = unique(x[, 2]),
        xvcount = tapply(x[, 3], x[, 2], sum),
        yvcount = tapply(y[, 3], y[, 2],sum))
    }
    if(smoothmeth != "none")
    {
      out$smoothmethod <- smoothmeth
      out$smoothout <- list(x = xsmooth, y = ysmooth)
    }
  }
  else out <- yx[, 1]
  return(out)
}
