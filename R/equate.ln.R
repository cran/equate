equate.ln <- function(x,y,type="linear",method="none",w=1,
  internal=TRUE,verbose=FALSE,...)
{
  xscale <- unique(x[,1])
  method <- match.arg(tolower(method),c("none","tucker","levine"))
  if(method=="none") 
  {
    mx <- mean(x)
    sdx <- sqrt(cov.freqtab(x))
    my <- mean(y)
    sdy <- sqrt(cov.freqtab(y))
    stats <- c(mx,sdx,my,sdy)
  }
  else
  {
    synthstats <- synthetic(x,y,w,method,internal)$s
    stats <- synthstats[c(5,11,6,12)]
  }

  type <- match.arg(tolower(type),c("mean","linear"))
  if(type=="mean") slope <- 1
  else slope <- stats[4]/stats[2]
  intercept <- stats[3]-slope*stats[1]
  yx <- slope*xscale+intercept

  out <- yx
  if(verbose)
  {
    out <- list(yx=yx)
    out$coefficients <- rbind(intercept,slope)
    if(method=="none") out$yx <- cbind(yx=yx,se=se.ln(x,y))
    else
    {
      out$synthstats <- synthstats[-(1:2),]
      out$anchortab <- cbind(scale=unique(x[,2]),
        xvcount=tapply(x[,3],x[,2],sum),
        yvcount=tapply(y[,3],y[,2],sum))
    }
  }
  return(out)
}