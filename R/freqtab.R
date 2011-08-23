freqtab <- function(xscale, x, vscale, v) {

  if(missing(v))
      ftab <- cbind(x = xscale,
        count = as.vector(table(factor(x, levels = xscale))))

  else {
    temptab <- table(factor(v, levels = vscale),
      factor(x, levels = xscale))
    ftab <- cbind(x = rep(xscale, each = length(vscale)),
      v = rep(vscale, length(xscale)), count = as.vector(temptab))
  }

  rownames(ftab) <- NULL
  class(ftab) <- "freqtab"

  return(ftab)
}
