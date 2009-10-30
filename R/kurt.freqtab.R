kurt.freqtab <- function(freqtab)
{
  sum(((freqtab[,1]-mean(freqtab))^4)*freqtab[,2])/
    (sum(freqtab[,2])-1)/var.freqtab(freqtab)^2
}
