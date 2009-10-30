equate <- function(x,y,scale,type,method="none",bootse=FALSE,
  internal=TRUE,Ky=max(scale),xv,yv,vscale,w,...)
{
  xtab <- freqtab(x,scale)
  ytab <- freqtab(y,scale)
  xcounts <- xtab[,2]
  ycounts <- ytab[,2]
  mx <- mean(xtab)
  sdx <- sqrt(var.freqtab(xtab))
  my <- mean(ytab)
  sdy <- sqrt(var.freqtab(ytab))

  type <- match.arg(tolower(type),c("mean","linear","equipercentile"))
  eqfun <- switch(type,equipercentile="equate.eq","equate.ln")
  eqfun <- match.fun(eqfun)
  eqout <- eqfun(x,y,scale,type=type,method=method,internal=internal,
    verbose=TRUE,Ky=Ky,xv=xv,yv=yv,vscale=vscale,w=w,...)

  yx <- freqtab(xcounts,eqout$yx[1:length(scale)])
  out <- list(type=type,method=match.arg(tolower(method),
    c("none","tucker","levine","frequency estimation")))
  out$design <- ifelse(method=="none","random groups","nonequivalent groups")
  out$stats <- cbind(mean=c(mx,my,mean(yx)),sd=c(sdx,sdy,
    sqrt(var.freqtab(yx))))
  rownames(out$stats) <- c("x","y","yx")
  out$freqtab <- cbind(scale=scale,xcounts=xcounts,fx=fx(xtab),
    ycounts=ycounts,fy=fx(ytab))
  out$concordance <- cbind(scale,eqout$yx)
  out <- c(out,eqout[-1])
  if(bootse)
  { 
    out$see <- se.boot(xtab,ytab,eqfun=eqfun,type=type,method=method,
      internal=internal,...)
  }
  class(out) <- "equate"
  return(out)
}