freqtab <- function(x,xscale,v,vscale,addclass=FALSE)
{
  if(missing(v))
  {
    if(addclass) freqtab <- cbind(xscale,x)
    else
      freqtab <- cbind(xscale,count=as.vector(table(factor(x,levels=xscale))))
  }
  else
  {
    temptab <- table(factor(v,levels=vscale),factor(x,levels=xscale))
    freqtab <- cbind(x=rep(xscale,each=length(vscale)),
      v=rep(vscale,length(xscale)),
      count=as.vector(temptab))
  }
  class(freqtab) <- "freqtab"
  return(freqtab)
}