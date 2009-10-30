se.boot <- function(x,y,xn=sum(x[,2]),yn=sum(y[,2]),reps=500,eqfun,...)
{
  scale <- x[,1]
  xprob <- x[,2]/sum(x[,2])
  yprob <- y[,2]/sum(y[,2])
  eqfun <- match.fun(eqfun)
  eqmat <- matrix(nrow=nrow(x),ncol=reps)
  for(i in 1:reps)
  {
    xtemp <- freqtab(sample(x[,1],xn,replace=TRUE,prob=xprob),x[,1])[,2]
    ytemp <- freqtab(sample(y[,1],yn,replace=TRUE,prob=yprob),y[,1])[,2]
    eqmat[,i] <- eqfun(xtemp,ytemp,scale,...)$yx
  }
  return(apply(eqmat,1,sd))
}