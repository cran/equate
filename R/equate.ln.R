equate.ln <- function(x,y,scale,type="linear",method="none",xv,yv,
  vscale,w=1,internal=TRUE,verbose=FALSE,...)
{
  xtab <- freqtab(x,scale)
  ytab <- freqtab(y,scale)
  mx <- mean(xtab)
  sdx <- sqrt(var.freqtab(xtab))
  my <- mean(ytab)
  sdy <- sqrt(var.freqtab(ytab))
  
  method <- match.arg(tolower(method),c("none","tucker","levine"))
  if(method=="none") stats <- c(mx,sdx,my,sdy)
  else
  {
    synthstats <- synthetic(x,xv,y,yv,w,method,internal)$s
    stats <- synthstats[c(5,11,6,12)]
  }

  type <- match.arg(tolower(type),c("mean","linear"))
  if(type=="mean") slope <- 1
  else slope <- stats[4]/stats[2]
  intercept <- stats[3]-slope*stats[1]
  yx <- slope*xtab[,1]+intercept

  out <- list(yx=cbind(yx))
  if(verbose)
  {
    out$coefficients <- rbind(intercept,slope)
    if(method=="none") out$yx <- cbind(yx=yx,se=se.ln(x,y,scale))
    else
    {
      out$synthstats <- synthstats[-(1:2),]
      out$anchortab <- cbind(scale=vscale,xvcount=freqtab(xv,vscale)[,2],
        yvcount=freqtab(yv,vscale)[,2])
    }
  }
  return(out)
}