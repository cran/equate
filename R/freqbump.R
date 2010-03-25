freqbump <- function(x, jmin = 10^-6, Kx = max(x[, 1])) {
  f <- x[, ncol(x)]
  if(sum(f) != 1) f <- f/sum(f)
  fbump = double()
  for(a in 1:length(f))
    fbump[a] <- (f[a] + jmin)/(1 + Kx * jmin)
  return(fbump)
}
