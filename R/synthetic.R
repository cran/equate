synthetic <- function(x, y, w = 1, method, internal = TRUE,
  lts = FALSE) {

  xscale <- unique(x[, 1])
  vscale <- unique(x[, 2])
  if(w == -1)
    w <- sum(x[, 3])/sum(x[, 3], y[, 3])
  yw <- 1 - w
  mx <- mean(x)
  varx <- cov.freqtab(x[, -2])
  mxv <- mean.freqtab(x[, -1])
  varxv <- cov.freqtab(x[, -1])
  my <- mean(y)
  vary <- cov.freqtab(y[, -2])
  myv <- mean.freqtab(y[, -1])
  varyv <- cov.freqtab(y[, -1])

  method <- match.arg(tolower(method), c("nominal weights",
    "tucker", "levine", "frequency estimation", "braun/holland"))
  if(method == "frequency estimation" |
    method == "braun/holland") {
    xj <- matrix(x[, 3]/sum(x[, 3]), ncol = length(vscale),
      nrow = length(xscale), byrow = TRUE)
    yj <- matrix(y[, 3]/sum(y[, 3]), ncol = length(vscale),
      nrow = length(xscale), byrow = TRUE)
    h1 <- apply(xj, 2, sum)
    h2 <- apply(yj, 2, sum)
    fx2 <- (xj %*% diag(1/h1)) %*% diag(h2)
    gy1 <- (yj %*% diag(1/h2)) %*% diag(h1)
    fs <- (w * apply(xj, 1, sum) + yw * apply(fx2, 1, sum)) *
      sum(x[, 3])
    gs <- (w * apply(gy1, 1, sum) + yw * apply(yj, 1, sum)) *
      sum(y[, 3])
    fstab <- as.freqtab(xscale, fs)
    gstab <- as.freqtab(xscale, gs)
    msx <- mean(fstab)
    msy <- mean(gstab)
    sdsx <- sqrt(cov.freqtab(fstab))
    sdsy <- sqrt(cov.freqtab(gstab))
  }
  else {
    covxv <- cov.freqtab(x)
    covyv <- cov.freqtab(y)

    if(method == "nominal weights")
      slope1 <- slope2 <- max(vscale)/max(xscale)
    if(method == "tucker") {
      slope1 <- covxv/varxv
      slope2 <- covyv/varyv
    }
    else if(method == "levine" & internal) {
      slope1 <- varx/covxv
      slope2 <- vary/covyv
    }
    else if(method == "levine" & !internal) {
      slope1 <- (varx + covxv)/(varxv + covxv)
      slope2 <- (vary + covyv)/(varyv + covyv)
    }
    if(!lts) {
      msx <- mx - (yw * slope1 * (mxv - myv))
      msy <- my + (w * slope2 * (mxv - myv))
      sdsx <- sqrt(varx -
        (yw * (slope1^2) * (varxv - varyv)) +
        (w * yw * (slope1^2) * (mxv - myv)^2))
      sdsy <- sqrt(vary +
        (w * (slope2^2) * (varxv - varyv)) +
        (w * yw * (slope2^2) * (mxv - myv)^2))
    }
  }
  if(lts)
    out <- list(gamma = c(slope1, slope2))
  else {
    out <- list(synthstats = rbind(xs = c(msx, sdsx),
      ys = c(msy, sdsy)))
    colnames(out$synthstats) <- c("mean", "sd")
    out$w <- w
  }
  if(method == "frequency estimation" |
    method == "braun/holland")
    out <- c(out, list(synthtab = cbind(xscale,
      xcount = fstab[, 2], ycount = gstab[, 2])))
  return(out)
}
