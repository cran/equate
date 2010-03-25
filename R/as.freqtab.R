as.freqtab <- function(xscale, x = NULL) {
  ftab <- cbind(xscale, x)
  class(ftab) <- "freqtab"
  if(ncol(ftab) == 3)
    colnames(ftab) <- c("x", "v", "count")
  else colnames(ftab) <- c("x", "count")
  rownames(ftab) <- NULL
  return(ftab)
}
