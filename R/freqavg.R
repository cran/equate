freqavg <- function(freqtab,jmin=1)
{
  freqtab <- cbind(freqtab,0,0,0,0,0)
  ks <- 1
  while(ks<=nrow(freqtab) & freqtab[ks,2]<jmin) ks <- ks+1
  freqtab[1:ks,3] <- 1
  freqtab[1:ks,4] <- ks
  
  ls=nrow(freqtab)
  while(ls>=0 & freqtab[ls,2]<jmin) ls <- ls-1
  freqtab[ls:nrow(freqtab),3] <- ls
  freqtab[ls:nrow(freqtab),4] <- nrow(freqtab)
  
  ss <- ks+1
  ts <- ls-1  
  for(j in ss:ts)
  {
    if(freqtab[j,2]<jmin)
    {
      ls <- j
      ks <- j
      while(ls>=1 & freqtab[ls,2]<jmin) ls <- ls-1
      while(ks<=nrow(freqtab) & freqtab[ks,2]<jmin) ks <- ks+1
      freqtab[ls:ks,3] <- ls
      freqtab[ls:ks,4] <- ks
    }
  }
  
  for(p in 1:(nrow(freqtab)-1))
  {
    if(freqtab[p,4]>0 & freqtab[p,4]==freqtab[p+1,3])
    {
      if(freqtab[p,3]>0) ls <- freqtab[p,3]
      if(freqtab[p+1,4]>0) ks <- freqtab[p+1,4]
      freqtab[ls:ks,3] <- ls
      freqtab[ls:ks,4] <- ks
    }
  }
  
  for(j in 1:nrow(freqtab))
  {
    ls <- freqtab[j,3]
      if(ls==0) ls <- j
    ks <- freqtab[j,4]
      if(ks==0) ks <- j
    sumit <- 0
    sumit <- sumit+sum(freqtab[ls:ks,2])
    for(i in ls:ks)
    {
      freqtab[i,5] <- sumit
      freqtab[i,6] <- freqtab[j,4]-freqtab[j,3]+1
      freqtab[i,7] <- freqtab[i,5]/freqtab[i,6]
    }
    j <- j+freqtab[j,4]-freqtab[j,3]
  }
  
  colnames(freqtab)[c(1,2,6,7)] <- c("score","count","b","acount")
  return(freqtab[,c(1,2,7)])
}