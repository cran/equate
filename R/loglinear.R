loglinear <- function(x,scorefun,degree,raw=TRUE,convergecrit=.0001,...)
{
  xscale <- unique(x[,1])
  if(missing(scorefun))
    scorefun <- poly(xscale,degree=degree,raw=raw)
  
  output <- list()
  counts <- x[,ncol(x)]
  colnames(scorefun) <- NULL
  b <- usb <- cbind(scorefun)
  ntot <- sum(counts)
  a <- .8*counts+.2*(ntot/length(counts)) 
  
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
  likold <- alpha*ntot+sum((t(counts)%*%b)%*%beta)
  a <- t(b)%*%diag(m)%*%b - t(b)%*%m%*%t(m)%*%b/ntot
  c <- t(b)%*%counts-t(b)%*%m
  delta <- solve(a,c)
  
  beta <- beta+delta
  alpha <- -log(sum(exp(b%*%beta)))
  lik <- alpha*ntot+sum((t(counts)%*%b)%*%beta)
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
    c <- t(b)%*%counts-t(b)%*%m
    delta <- solve(a,c)
    beta <- beta+delta
    alpha <- -log(sum(exp(b%*%beta)))
    m <- as.vector(ntot*exp(alpha+b%*%beta))
  
    ebc2 <- matrix(0,nrow=ncol(b),ncol=3)
    for(i in 1:ncol(b))
    {
      ebc2[i,1] <- t(b[,i])%*%counts
      ebc2[i,2] <- t(b[,i])%*%m
      if(ebc2[i,1]==0) ebc2[i,3] <- abs(ebc2[i,2])
      if(abs(ebc2[i,1])>0) ebc2[i,3] <- abs((ebc2[i,1]-ebc2[i,2])/ebc2[i,1])
    }
    testerc2 <- max(ebc2[,3])
    if(testerc2>convergecrit) crit2 <- 0 
      else crit2 <- 1
    lik <- alpha*ntot+sum((t(counts)%*%b)%*%beta)
  }
  
  liksum1 <- counts/m
  liksum1[liksum1==0] <- 1
  liksum2 <- log(liksum1)
  likchisquare <- 2*(t(counts)%*%liksum2)
  aic <- likchisquare+2*(ncol(b)+1)
  caic <- likchisquare+(1+log(ntot))*(ncol(b)+1)
  pearsonchisquare <- sum((counts-m)^2/m)
  freemantukeychisquare <- sum((sqrt(counts)+sqrt(counts+1)-sqrt(4*m+1))^2)
  df <- length(counts)-ncol(b)-1
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
  
  ftres <- sqrt(counts)+sqrt(counts+1)-sqrt(4*m+1)
  output$fitted.values <- m
  output$residuals <- ftres
  output$cmatrix  <- C
  colnames(output$cmatrix) <- colnames(scorefun) #paste("C.",1:ncol(b),sep="")
  output$scorefun <- scorefun
  return(output)
}