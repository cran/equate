fx <- function(x) {
  colnames(x) = NULL
  x[, 2] <- x[, 2]/sum(x[, 2])
  f <- vector()
  for(i in 1:nrow(x)) f[i] <- sum(x[1:i, 2])
  return(f)
}
