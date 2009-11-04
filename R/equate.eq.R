equate.eq <- function(x,y,method="none",Ky=max(y[,1]),w=1,smooth="none",
  jmin,xscorefun,yscorefun,verbose=FALSE,...)
{
  nc <- ncol(x)
  xscale <- unique(x[,1])
  xtab <- x
  ytab <- y
  smooth <- match.arg(tolower(smooth),c("none","bump","average","loglin"))
  if(smooth=="bump")
  {
    if(missing(jmin)) jmin <- 10^-6
    xsmooth <- freqbump(x[,nc],jmin,Ky)*sum(x[,nc])
    ysmooth <- freqbump(y[,nc],jmin,Ky)*sum(y[,nc])
    xtab[,nc] <- xsmooth
    ytab[,nc] <- ysmooth
  }
  else if(smooth=="average" & nc==2)
  {
    if(missing(jmin)) jmin <- 1
    xsmooth <- freqavg(x,jmin) 
    ysmooth <- freqavg(y,jmin)
    xtab[,2] <- xsmooth
    ytab[,2] <- ysmooth
  }
  else if(smooth=="loglin")
  {
    xsmooth <- loglinear(x,xscorefun,...)
    ysmooth <- loglinear(y,yscorefun,...)
    xtab[,nc] <- xsmooth$fitted
    ytab[,nc] <- ysmooth$fitted
  }
  method <- match.arg(tolower(method),c("none","frequency"))
  if(method=="frequency" | nc==3)
  {
    stabs <- synthetic(xtab,ytab,w,method)
    xtab <- cbind(xscale,stabs$synthtab[,2])
    ytab <- cbind(xscale,stabs$synthtab[,3])
  }

  xtab <- cbind(xtab,px=px(xtab),0)
  ytab <- cbind(ytab,fx=fx(ytab))
  sn <- nrow(xtab)
  se <- vector(length=sn)
  numK <- is.numeric(Ky)
  xn <- sum(xtab[,2])
  yn <- sum(ytab[,2])
  Ly <- min(xscale)
  xone <- which(xtab[,3]==1)
  xless <- sum(xtab[,3]!=1)

  if(numK) xtab[xone,4] <- Ky+.5
  else xtab[xone,4] <- xtab[xone,1]
  for(j in 1:xless)
  {
    yu <- 1
    while(ytab[yu,3]<=xtab[j,3]) yu <- yu+1
    yu2 <- ifelse(yu==1, 0, ytab[yu-1,3])
    g0 <- ytab[yu,3]-yu2
    xtab[j,4] <- ytab[yu,1]-.5+((xtab[j,3]-yu2)/g0)
    if(g0>0) se[j] <- se.eq(xtab[j,3],g0,yu2,xn,yn)
  }
  if(any(ytab[,2]==0))
  {
    for(i in 1:xless)
    {
      if(xtab[i,3]==0) xtab[i,4] <- ifelse(numK,Ly-.5,xtab[i,1])
      else
      {
        yl <- sn
        while(ytab[yl,3]>=xtab[i,3]) yl <- yl-1
        yl2 <- ytab[yl+1,3]
        xtab.temp <- ytab[yl,1]+.5+((xtab[i,3]-ytab[yl,3])/(yl2-ytab[yl,3]))
        xtab[i,4] <- mean(c(xtab[i,4],xtab.temp))
      }
    }
  }
  
  out <- xtab[,4]
  if(verbose)
  {
    out <- list(yx=xtab[,4])
    if(smooth!="none")
    {
      out$smoothmethod <- smooth 
      out$smoothout <- list(x=xsmooth,y=ysmooth)
    }
    if(method=="frequency")
    {
      out <- c(out,stabs)
      out$anchortab <- cbind(scale=unique(x[,2]),
        xvcount=tapply(x[,3],x[,2],sum),
        yvcount=tapply(y[,3],y[,2],sum))
    }
    else out$yx <- cbind(yx=xtab[,4],se)
  }
  return(out)
}