equate.ca <- function(x, y, method = NA, lowp, midp = "mean",
  highp, verbose = FALSE, ...) {

  xscale <- unique(x[, 1])
  yscale <- unique(y[, 1])

  x1 <- ifelse(missing(lowp), min(xscale), lowp[1])
  y1 <- ifelse(missing(lowp), min(yscale), lowp[length(lowp)])
  x3 <- ifelse(missing(highp), max(xscale), highp[1])
  y3 <- ifelse(missing(highp), max(yscale), highp[length(highp)])
  index <- xscale >= x1 & xscale <= x3
  slope <- (y3 - y1)/(x3 - x1)
  intercept <- y1 - slope * x1
  lin <- cyx <- slope * xscale + intercept

  x2 <- mean(x)
  method <- match.arg(tolower(method),
    c(NA, "nominal weights", "tucker", "levine", "chained",
      "braun/holland"))
  if(is.na(method))
    y2 <- mean(y)
  else if(method == "chained") {
    sdy <- sd.freqtab(y[, -2])
    sdyv <- sd.freqtab(y[, -1])
    mxv <- mean.freqtab(x[, -1])
    myv <- mean.freqtab(y[, -1])
    slope2 <- ifelse(midp == "mean", 1, sdy/sdyv)
    y2 <- mean(y) + slope2 * (mxv - myv)
  }
  else {
    lyx <- equate(x, y, type = "mean", method = method, ...)
    y2 <- equate(x2, lyx)
  }
  y2star <- y2 - x2

  xcent <- (x3^2 - x1^2)/(2 * (x3 - x1))
  ycent <- ((x1^2) * (x3 - x2) - (x2^2 + y2star^2) *
    (x3 - x1) + (x3^2) * (x2 - x1))/(2 * (y2star * (x1 - x3)))
  r <- sqrt((xcent - x1)^2 + ycent^2)

  if(y2star < 0)
    cyx[index] <- lin[index] + ycent -
      sqrt((r^2) - (lin[index] - xcent)^2)
  else
    cyx[index] <- lin[index] + ycent +
      sqrt((r^2) - (lin[index] - xcent)^2)
  out <- cyx
  if(verbose) {
    out <- list(yx = cyx)
    out$coefficients <- rbind(intercept, slope, xcenter = xcent,
      ycenter = ycent, r)[,1]
    out$points <- rbind(x = c(x1, x2, x3), y = c(y1, y2, y3))
    colnames(out$points) <- c("low", "middle", "high")
    if(!is.na(method)) {
      out$anchorstats <- rbind(descript(x[, -1]),
        descript(y[, -1]))
      rownames(out$anchorstats) <- c("xv", "yv")
      out$anchortab <- cbind(scale = unique(x[, 2]),
        xvcount = tapply(x[, 3], x[, 2], sum),
        yvcount = tapply(y[, 3], y[, 2], sum))
    }
  }
  return(out)
}
