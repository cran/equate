px <- function(freqtab)
{
  colnames(freqtab)=NULL
  freqtab[,2] <- freqtab[,2]/sum(freqtab[,2])
  px <- .5*freqtab[1,2]
  for(i in 2:nrow(freqtab)) px[i] <- sum(freqtab[1:i-1,2])+.5*freqtab[i,2]
  return(px)
}