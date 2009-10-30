freqtab <- function(x,scale)
{
  if(length(x)==length(scale)) freqtab <- cbind(scale,x)
  else freqtab <- cbind(scale,count=as.vector(table(factor(x,levels=scale))))
  class(freqtab) <- "freqtab"
  return(freqtab)
}