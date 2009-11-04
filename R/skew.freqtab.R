skew.freqtab <- function(x)
{
  sum(((x[,1]-mean(x))^3)*x[,2])/(sum(x[,2])-1)/cov.freqtab(x)^1.5
}
