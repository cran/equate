fx <- function(x)
{
  colnames(x)=NULL
  x[,2] <- x[,2]/sum(x[,2])
  fx <- vector()
  for(i in 1:nrow(x)) fx[i] <- sum(x[1:i,2])
  return(fx)
}