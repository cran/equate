var.freqtab <- function(freqtab)
{
  sum(((freqtab[,1]-mean(freqtab))^2)*freqtab[,2])/(sum(freqtab[,2])-1)
}