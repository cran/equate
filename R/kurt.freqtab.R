kurt.freqtab <- function(x)
{
  sum(((x[,1]-mean(x))^4)*x[,2])/(sum(x[,2])-1)/cov.freqtab(x)^2
}
