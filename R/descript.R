descript <- function(x) {
  return(cbind(mean = mean.freqtab(x),
    sd = sd.freqtab(x),
    skew = skew.freqtab(x),
    kurt = kurt.freqtab(x),
    n = sum(x[, ncol(x)])))
}
