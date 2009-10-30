equate.eq <- function(x,y,scale,method="none",Ky=max(scale),xv,yv,
  vscale,w,smooth="none",jmin,xscorefun,yscorefun,verbose=FALSE,...)
{
  xtab <- freqtab(x,scale)
  ytab <- freqtab(y,scale)

  smooth <- match.arg(tolower(smooth),c("none","bump","average","loglin"))
  if(smooth=="bump")
  {
    if(missing(jmin)) jmin <- 10^-6
    xsmooth <- freqbump(xtab[,2],jmin,Ky)*sum(xtab[,2])
    ysmooth <- freqbump(ytab[,2],jmin,Ky)*sum(ytab[,2])
    xtab[,2] <- xsmooth
    ytab[,2] <- ysmooth
  }
  else if(smooth=="average")
  {
    if(missing(jmin)) jmin <- 1
    xsmooth <- freqavg(xtab,jmin) 
    ysmooth <- freqavg(ytab,jmin)
    xtab <- xsmooth[,c(1,3)]
    ytab <- ysmooth[,c(1,3)]
  }
  else if(smooth=="loglin")
  {
    xsmooth <- loglinear(x,scale,xscorefun,...)
    ysmooth <- loglinear(y,scale,yscorefun,...)
    xtab[,2] <- xsmooth$fitc[,"smoothedcounts"]
    ytab[,2] <- ysmooth$fitc[,"smoothedcounts"]
  }
  
  method <- match.arg(tolower(method),c("none","frequency"))
  if(method=="frequency")
  {
    stabs <- synthetic(x,xv,y,yv,w,method,scale=scale,
      vscale=vscale)
    xtab[,2] <- stabs$synthtab[,2]
    ytab[,2] <- stabs$synthtab[,3]
  }

  xtab <- cbind(xtab,px=px(xtab),0)
  ytab <- cbind(ytab,fx=fx(ytab))
  sn <- nrow(xtab)
  se <- vector(length=sn)
  numK <- is.numeric(Ky)
  xn <- sum(xtab[,2])
  yn <- sum(ytab[,2])
  Ly <- min(scale)
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
  
  out <- list(yx=xtab[,4])
  if(verbose)
  {
    if(smooth!="none")
    {
      out$smoothmethod <- smooth 
      out$smoothout <- list(xsmooth=xsmooth,ysmooth=ysmooth)
    }
    if(method=="frequency")
    {
      out <- c(out,stabs)
      out$anchortab <- cbind(scale=vscale,xvcount=freqtab(xv,vscale)[,2],
        yvcount=freqtab(yv,vscale)[,2])
    }
    else out$yx <- cbind(yx=xtab[,4],se)
  }
  return(out)
}