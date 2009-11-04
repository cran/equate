se.boot <- function(x,y,xn=sum(x[,ncol(x)]),yn=sum(y[,ncol(y)]),
  reps=500,eqfun,...)
{
  nc <- ncol(x)
  xscale <- unique(x[,1])
  xprob <- x[,nc]/sum(x[,nc])
  yprob <- y[,nc]/sum(y[,nc])
  eqfun <- match.fun(eqfun)
  eqmat <- matrix(nrow=length(xscale),ncol=reps)
  if(nc==3)
  {
    index <- 1:nrow(x)
    vscale <- unique(x[,2])
    for(i in 1:reps)
    {
      xi <- sample(index,xn,replace=TRUE,prob=xprob)
      xtemp <- freqtab(x[index,1],xscale,x[index,2],vscale)
      yi <- sample(index,yn,replace=TRUE,prob=yprob)
      ytemp <- freqtab(y[index,1],xscale,y[index,2],vscale)
      eqmat[,i] <- eqfun(xtemp,ytemp,...)
    }
  }
  else
  {
    for(i in 1:reps)
    {
      xtemp <- freqtab(sample(x[,1],xn,replace=TRUE,prob=xprob),xscale)
      ytemp <- freqtab(sample(y[,1],yn,replace=TRUE,prob=yprob),xscale)
      eqmat[,i] <- eqfun(xtemp,ytemp,...)
    }
  }
  return(apply(eqmat,1,sd))
}