equipercentile <- function(x, y, Ky = max(y[, 1])) {
  if(Ky == "maxobs")
    Ky <- max.freqtab(y)
  yscale <- unique(y[, 1])
  yn <- sum(y[, 2])
  if(class(x) != "freqtab") {
    prank <- x
    xscale <- yscale
    xn <- 0
  }
  else {
    prank <- px(x)
    xscale <- unique(x[, 1])
    xn <- sum(x[, 2])
    se <- vector(length = length(prank))
  }
  yx <- vector()
  fy <- fx(y)
  sn <- nrow(y)
  Ly <- min(yscale)
  xnone <- prank == 0
  xone <- prank == 1
  xbot <- sum(xnone) + 1
  xtop <- sum(!xone)
  xyone <- which(xscale[xone] > (Ky + .5)) + xtop

  yx[xnone] <- Ly - .5
  yx[xone] <- Ky + .5
  yx[xyone] <- xscale[xyone]
  for(j in xbot:xtop) {
    yu <- 1
    while(fy[yu] <= prank[j]) yu <- yu + 1
    yu2 <- ifelse(yu == 1, 0, fy[yu - 1])
    g0 <- fy[yu] - yu2
    yx[j] <- y[yu, 1] - .5 + ((prank[j] - yu2)/g0)
    if(g0 > 0 & xn)
      se[j] <- se.eq(prank[j], g0, yu2, xn, yn)
  }
  if(any(y[, 2] == 0)) {
    xbot <- xbot + sum(prank[!xnone] <= min(fy))
    for(i in xbot:xtop) {
      yl <- sn
      while(fy[yl] >= prank[i]) yl <- yl - 1
      yl2 <- fy[yl + 1]
      yxtemp <- y[yl, 1] + .5 + ((prank[i] - fy[yl])/
        (yl2 - fy[yl]))
      yx[i] <- mean(c(yx[i], yxtemp))
    }
  }

  out <- list(yx = yx)
  if(xn)
    out$se <- se
  return(out)
}
