loglinear <- function(x,scale,scorefun,degree,raw=TRUE,convergecrit=.0001,...)
{
  if(missing(scorefun))
    scorefun <- poly(scale,degree=degree,raw=raw)
  
  freqtab <- freqtab(x,scale)
  output <- list()
  n <- freqtab[,2]
  colnames(scorefun) <- NULL
  b <- usb <- cbind(scorefun)
  ntot <- sum(n)
  a <- .8*n+.2*(ntot/length(n)) 
  
  bssq <- rep(0,ncol(b))
  for(i in 1:ncol(b))
  {
    bssq[i] <- sum(b[,i]^2)
    if(sum(b[,i])!=nrow(b))
    {
      b[,i] <- b[,i]-sum(b[,i])/nrow(b)
      b[,i] <- b[,i]/sqrt(bssq[i])
    }
  }
  
  ibeta1 <- t(b)%*%diag(a)%*%b - t(b)%*%a%*%t(a)%*%b/ntot
  ibeta2 <- t(b)%*%diag(a)%*%log(a) - t(b)%*%a%*%t(a)%*%log(a)/ntot
  beta <- as.matrix(solve(ibeta1,ibeta2))
  
  alpha <- -log(sum(exp(b%*%beta)))
  m <- as.vector(ntot*exp(alpha+b%*%beta))
  likold <- alpha*ntot+sum((t(n)%*%b)%*%beta)
  a <- t(b)%*%diag(m)%*%b - t(b)%*%m%*%t(m)%*%b/ntot
  c <- t(b)%*%n-t(b)%*%m
  delta <- solve(a,c)
  
  beta <- beta+delta
  alpha <- -log(sum(exp(b%*%beta)))
  lik <- alpha*ntot+sum((t(n)%*%b)%*%beta)
  iter <- 1
  crit1 <- 0
  crit2 <- 0
  
  while(crit1+crit2<2)
  {
    if (abs((lik-likold)/likold)<convergecrit) crit1 <- 1
    likchang <- abs((lik-likold)/likold)
    likold <- lik
    iter <- iter+1
    a <- t(b)%*%diag(m)%*%b - t(b)%*%m%*%t(m)%*%b/ntot
    c <- t(b)%*%n-t(b)%*%m
    delta <- solve(a,c)
    beta <- beta+delta
    alpha <- -log(sum(exp(b%*%beta)))
    m <- as.vector(ntot*exp(alpha+b%*%beta))
  
    ebc2 <- matrix(0,nrow=ncol(b),ncol=3)
    for(i in 1:ncol(b))
    {
      ebc2[i,1] <- t(b[,i])%*%n
      ebc2[i,2] <- t(b[,i])%*%m
      if(ebc2[i,1]==0) ebc2[i,3] <- abs(ebc2[i,2])
      if(abs(ebc2[i,1])>0) ebc2[i,3] <- abs((ebc2[i,1]-ebc2[i,2])/ebc2[i,1])
    }
    testerc2 <- max(ebc2[,3])
    if(testerc2>convergecrit) crit2 <- 0 
      else crit2 <- 1
    lik <- alpha*ntot+sum((t(n)%*%b)%*%beta)
  }
  
  liksum1 <- n/m
  liksum1[liksum1==0] <- 1
  liksum2 <- log(liksum1)
  likchisquare <- 2*(t(n)%*%liksum2)
  aic <- likchisquare+2*(ncol(b)+1)
  caic <- likchisquare+(1+log(ntot))*(ncol(b)+1)
  pearsonchisquare <- sum((n-m)^2/m)
  freemantukeychisquare <- sum((sqrt(n)+sqrt(n+1)-sqrt(4*m+1))^2)
  df <- length(n)-ncol(b)-1
  output$modelfit <- rbind(likchisquare,pearsonchisquare,freemantukeychisquare,aic,caic,df)
    rownames(output$modelfit) <- c('Likelihood Ratio Chi-square','Pearson Chi-square',
      'Freeman-Tukey Chi-square','AIC','CAIC','Degrees of Freedom')
  
  a <- t(usb)%*%diag(m)%*%usb - t(usb)%*%m%*%t(m)%*%usb/ntot
  tscorefunctions <- t(scorefun)
  beta2 <- beta/sqrt(bssq)
  stdbeta <- tryCatch(sqrt(diag(solve(a,...))),error = function(x) NA)
  output$rawbetas <- cbind(beta2,stdbeta)
    colnames(output$rawbetas) <- c("beta","se")
  output$alpha <- -log(sum(exp(usb%*%beta2)))
  output$iterations <- iter
  
  mf <- m/sum(m)
  sqm <- sqrt(mf)
  D <- diag(sqm)
  ones <- matrix(1,nrow=nrow(b),ncol=1)
  
  qrsum <- ones%*%t(ones)%*%diag(mf)%*%b
  QR <- sqrt(mf)*(b-qrsum)
  
  Q <- qr(QR)$qr[,1:ncol(b)]
  Nt <- 1/sqrt(ntot)
  C <- Nt*D%*%Q
  
  p <- m/sum(m)
  ftres <- sqrt(n)+sqrt(n+1)-sqrt(4*m+1)
  output$fitc <- cbind(m,p,ftres,C)
  colnames(output$fitc) <- c("smoothedcounts","smoothedprobs",
    "ftresiduals",paste("C.",1:ncol(b),sep=""))
  output$scorefun <- scorefun
  return(output)
}