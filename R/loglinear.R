loglinear <- function(x, scorefun, degree, raw = TRUE,
  convergecrit = .0001, verbose = TRUE, ...) {

  a.fn <- function(ab, am, antot){
    atemp <- matrix(nrow = ncol(ab), ncol = ncol(ab))
    for(i in 1:ncol(ab)){
      for(j in 1:ncol(ab)){
        ba1 <- sum(ab[, i] * am)
        ba2 <- sum(ab[, j] * am)
        ba3 <- sum(ab[, i] * ab[, j] * am)
        atemp[i, j] <- ba3 - (ba1 * ba2/antot)
      }
    }
    return(atemp)
  }

  nc <- ncol(x)
  if(missing(scorefun)) {
    scorefun <- vector()
    for(i in 1:(nc - 1))
      scorefun <- cbind(scorefun, poly(x[, i],
        degree = degree, raw = raw))
  }

  output <- list()
  counts <- x[, nc]
  colnames(scorefun) <- NULL
  b <- usb <- cbind(scorefun)
  ntot <- sum(counts)
  a <- .8 * counts + .2 * (ntot/length(counts))

  bssq <- rep(0, ncol(b))
  for(i in 1:ncol(b)) {
    bssq[i] <- sum(b[, i]^2)
    if(sum(b[, i]) != nrow(b)) {
      b[, i] <- b[, i] - sum(b[, i])/nrow(b)
      b[, i] <- b[, i]/sqrt(bssq[i])
    }
  }

  iB1 <- matrix(nrow = ncol(b), ncol = ncol(b))
  for(k in 1:ncol(b)){
    for(j in 1:ncol(b)){
      ba1 <- ba2 <- ba3 <- 0
      for(i in 1:nrow(b)){
        ba1 <- ba1 + b[i, k] * a[i]
        ba2 <- ba2 + b[i, j] * a[i]
        ba3 <- ba3 + b[i, k] * b[i, j] * a[i]
      }
      iB1[k, j] <- ba3 - (ba1 * ba2/ntot)
    }
  }

  iB2 <- vector(length = ncol(b))
  for(i in 1:ncol(b)){
    ba1 <- sum(b[, i] * a)
    ba2 <- sum(log(a) * a)
    ba3 <- sum(b[, i] * log(a) * a)
    iB2[i] <- ba3 - (ba1 * ba2/ntot)
  }

  B <- as.matrix(solve(iB1, iB2))

  alpha <- -log(sum(exp(b %*% B)))
  m <- as.vector(ntot * exp(alpha + b %*% B))
  likold <- alpha * ntot + sum((t(counts) %*% b) %*% B)
  a <- a.fn(b, m, ntot)
  c <- t(b) %*% counts - t(b) %*% m
  delta <- solve(a, c)

  B <- B + delta
  alpha <- -log(sum(exp(b %*% B)))
  lik <- alpha * ntot + sum((t(counts) %*% b) %*% B)
  iter <- 1
  crit1 <- 0
  crit2 <- 0

  while(crit1 + crit2 < 2) {
    if (abs((lik - likold)/likold) < convergecrit)
      crit1 <- 1
    likchang <- abs((lik - likold)/likold)
    likold <- lik
    iter <- iter + 1
    a <- a.fn(b, m, ntot)
    c <- t(b) %*% counts - t(b) %*% m
    delta <- solve(a, c)
    B <- B + delta
    alpha <- -log(sum(exp(b %*% B)))
    m <- as.vector(ntot * exp(alpha + b %*% B))

    ebc2 <- matrix(0, nrow = ncol(b), ncol = 3)
    for(i in 1:ncol(b)) {
      ebc2[i,1] <- t(b[, i]) %*% counts
      ebc2[i,2] <- t(b[, i]) %*% m
      if(ebc2[i, 1] == 0)
        ebc2[i, 3] <- abs(ebc2[i, 2])
      if(abs(ebc2[i, 1]) > 0)
        ebc2[i, 3] <- abs((ebc2[i, 1] - ebc2[i, 2])/ebc2[i, 1])
    }
    testerc2 <- max(ebc2[, 3])
    if(testerc2 > convergecrit)
      crit2 <- 0
    else crit2 <- 1
    lik <- alpha * ntot + sum((t(counts) %*% b) %*% B)
  }

  if(verbose) {
    liksum1 <- counts/m
    liksum1[liksum1 == 0] <- 1
    liksum2 <- log(liksum1)
    likchisquare <- 2 * (t(counts) %*% liksum2)
    aic <- likchisquare + 2 * (ncol(b) + 1)
    caic <- likchisquare + (1 + log(ntot)) * (ncol(b) + 1)
    pearsonchisquare <- sum((counts - m)^2/m)
    freemantukeychisquare <-
      sum((sqrt(counts) + sqrt(counts + 1) - sqrt(4 * m + 1))^2)
    dfree <- length(counts) - ncol(b) - 1
    output$modelfit <- rbind(likchisquare, pearsonchisquare,
      freemantukeychisquare, aic, caic, dfree)
    rownames(output$modelfit) <- c("Likelihood Ratio Chi-square",
      "Pearson Chi-square", "Freeman-Tukey Chi-square",
      "AIC", "CAIC", "Degrees of Freedom")

    a <- a.fn(usb, m, ntot)
    B2 <- B/sqrt(bssq)
    stdbeta <- tryCatch(sqrt(diag(solve(a, ...))),
      error = function(x) NA)
    output$rawbetas <- cbind(B2, stdbeta)
    colnames(output$rawbetas) <- c("beta", "se")
    output$alpha <- -log(sum(exp(usb %*% B2)))
    output$iterations <- iter

    mf <- m/sum(m)
    sqm <- sqrt(mf)
    ones <- matrix(1, nrow = nrow(b), ncol = 1)

    QR <- matrix(nrow = nrow(b), ncol = ncol(b))
    for(i in 1:ncol(b)){
      qrsum <- 0
      for(j in 1:nrow(b))
        qrsum <- qrsum + b[j, i] * mf[j]
      for(k in 1:nrow(b))
        QR[k, i] <- sqrt(mf[k]) * (b[k, i] - qrsum)
    }
    Q <- qr(QR)$qr[, 1:ncol(b)]
    Nt <- 1/sqrt(ntot)
    Cmat <- matrix(nrow = nrow(b), ncol = ncol(b))
    for(i in 1:nrow(b))
      Cmat[i,] <- Nt * sqm[i] * Q[i,]

    ftres <- sqrt(counts) + sqrt(counts + 1) - sqrt(4 * m + 1)
    output$fitted.values <- m
    output$residuals <- ftres
    output$cmatrix  <- Cmat
    colnames(output$cmatrix) <- colnames(scorefun)
    output$scorefun <- scorefun
  }
  else output$fitted.values <- m
  return(output)
}
