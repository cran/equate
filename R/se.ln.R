se.ln <- function(x,y,scale)
{
  x <- freqtab(x,scale)
  y <- freqtab(y,scale)
  nx <- sum(x[,2])
  ny <- sum(y[,2])
  mx <- mean(x)
  sdx <- sqrt(var.freqtab(x))
  vary <- var.freqtab(y)
  skewterm <- skew.freqtab(x)/nx+skew.freqtab(y)/ny
  kurtterm <- (kurt.freqtab(x)-1)/(4*nx)+(kurt.freqtab(y)-1)/(4*ny)
  se <- vector()
  for(i in 1:nrow(x))
  {
    xterm <- (scale[i]-mx)/sdx
    se[i] <- vary*(1/nx+1/ny+skewterm*xterm+kurtterm*xterm^2)
  }
  return(se)
}
