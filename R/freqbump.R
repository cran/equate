freqbump <- function(x,jmin=10^-6,Kx=length(x))
{
  if(!is.numeric(Kx)) Kx <- length(x)
  if(sum(x)!=1) x <- x/sum(x)
  fbump=double()
  for(a in 1:length(x)) fbump[a] <- (x[a]+jmin)/(1+Kx*jmin)
  return(fbump)
}