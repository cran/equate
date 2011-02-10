fx <- function(x) {
  if(!is.null(dim(x)))
    x <- x[, ncol(x)]
  x <- x/sum(x)
  f <- vector()
  for(i in 1:length(x)) f[i] <- sum(x[1:i])
  return(f)
}
