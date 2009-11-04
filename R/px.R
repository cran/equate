px <- function(x)
{
  colnames(x)=NULL
  x[,2] <- x[,2]/sum(x[,2])
  px <- .5*x[1,2]
  for(i in 2:nrow(x)) px[i] <- sum(x[1:i-1,2])+.5*x[i,2]
  return(px)
}