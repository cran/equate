fx <- function(freqtab)
{
  colnames(freqtab)=NULL
  freqtab[,2] <- freqtab[,2]/sum(freqtab[,2])
  fx <- vector()
  for(i in 1:nrow(freqtab)) fx[i] <- sum(freqtab[1:i,2])
  return(fx)
}