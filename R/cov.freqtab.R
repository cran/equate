cov.freqtab <- function(x)
{
  nc <- ncol(x)
  i <- x[,nc]>0
  xc <- x[i,1]-mean(x)
  vc <- x[i,nc-1]-mean.freqtab(x[,(nc-1):nc])
  return(sum(xc*vc*x[i,nc])/(sum(x[,nc])-1))
}