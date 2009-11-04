freqavg <- function(x,jmin=1)
{
  x <- cbind(x,0,0,0,0,0)
  ks <- 1
  while(ks<=nrow(x) & x[ks,2]<jmin) ks <- ks+1
  x[1:ks,3] <- 1
  x[1:ks,4] <- ks
  
  ls=nrow(x)
  while(ls>=0 & x[ls,2]<jmin) ls <- ls-1
  x[ls:nrow(x),3] <- ls
  x[ls:nrow(x),4] <- nrow(x)
  
  ss <- ks+1
  ts <- ls-1  
  for(j in ss:ts)
  {
    if(x[j,2]<jmin)
    {
      ls <- j
      ks <- j
      while(ls>=1 & x[ls,2]<jmin) ls <- ls-1
      while(ks<=nrow(x) & x[ks,2]<jmin) ks <- ks+1
      x[ls:ks,3] <- ls
      x[ls:ks,4] <- ks
    }
  }
  
  for(p in 1:(nrow(x)-1))
  {
    if(x[p,4]>0 & x[p,4]==x[p+1,3])
    {
      if(x[p,3]>0) ls <- x[p,3]
      if(x[p+1,4]>0) ks <- x[p+1,4]
      x[ls:ks,3] <- ls
      x[ls:ks,4] <- ks
    }
  }
  
  for(j in 1:nrow(x))
  {
    ls <- x[j,3]
      if(ls==0) ls <- j
    ks <- x[j,4]
      if(ks==0) ks <- j
    sumit <- 0
    sumit <- sumit+sum(x[ls:ks,2])
    for(i in ls:ks)
    {
      x[i,5] <- sumit
      x[i,6] <- x[j,4]-x[j,3]+1
      x[i,7] <- x[i,5]/x[i,6]
    }
    j <- j+x[j,4]-x[j,3]
  }
  
  colnames(x)[c(1,2,6,7)] <- c("score","count","b","acount")
  return(x[,7])
}