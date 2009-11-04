equate <- function(x,y,type,method="none",bootse=FALSE,
  internal=TRUE,Ky=max(y[,1]),w=1,...)
{
  nc <- ncol(x)
  xscale <- unique(x[,1])
  xcounts <- tapply(x[,nc],x[,1],sum)
  ycounts <- tapply(y[,nc],y[,1],sum)
  xtab <- freqtab(xcounts,xscale,addclass=TRUE)
  ytab <- freqtab(ycounts,xscale,addclass=TRUE)
  mx <- mean(xtab)
  sdx <- sqrt(cov.freqtab(xtab))
  my <- mean(ytab)
  sdy <- sqrt(cov.freqtab(ytab))

  type <- match.arg(tolower(type),c("mean","linear","equipercentile"))
  eqfun <- switch(type,equipercentile="equate.eq","equate.ln")
  eqfun <- match.fun(eqfun)
  eqout <- eqfun(x,y,type=type,method=method,internal=internal,
    verbose=TRUE,Ky=Ky,w=w,...)

  yx <- freqtab(xcounts,eqout$yx[1:length(xscale)],addclass=TRUE)
  out <- list(type=type,method=match.arg(tolower(method),
    c("none","tucker","levine","frequency estimation")))
  out$design <- ifelse(method=="none","random groups","nonequivalent groups")
  out$stats <- cbind(mean=c(mx,my,mean(yx)),sd=c(sdx,sdy,
    sqrt(cov.freqtab(yx))))
  rownames(out$stats) <- c("x","y","yx")
  out$freqtab <- cbind(scale=xscale,xcounts=xcounts,fx=fx(xtab),
    ycounts=ycounts,fy=fx(ytab))
  out$concordance <- cbind(scale=xscale,yx=eqout$yx)
  out <- c(out,eqout[-1])
  if(bootse)
  { 
    out$see <- se.boot(x,y,eqfun=eqfun,type=type,method=method,
      internal=internal,...)
  }
  class(out) <- "equate"
  return(out)
}