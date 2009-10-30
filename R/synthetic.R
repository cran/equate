synthetic <- function(x,xv,y,yv,w,method,internal=TRUE,scale,vscale)
{
  yw <- 1-w
  mx <- mean(x)
  varx <- var(x)
  mxv <- mean(xv)
  varxv <- var(xv)
  my <- mean(y)
  vary <- var(y)
  myv <- mean(yv)
  varyv <- var(yv)

  method <- match.arg(tolower(method),c("tucker","levine","frequency"))  
  if(method=="frequency")
  {
    xj <- table(factor(x,levels=scale),factor(xv,levels=vscale))/length(x)
    yj <- table(factor(y,levels=scale),factor(yv,levels=vscale))/length(y)
    h1 <- apply(xj,2,sum)
    h2 <- apply(yj,2,sum)
    fx2 <- (xj%*%diag(1/h1))%*%diag(h2)
    gy1 <- (yj%*%diag(1/h2))%*%diag(h1)
    fs <- (w*apply(xj,1,sum)+yw*apply(fx2,1,sum))*length(x)
    gs <- (w*apply(gy1,1,sum)+yw*apply(yj,1,sum))*length(y)
    fstab <- freqtab(fs,scale)
    gstab <- freqtab(gs,scale)
    msx <- mean(fstab)
    msy <- mean(gstab)
    sdsx <- sqrt(var.freqtab(fstab))
    sdsy <- sqrt(var.freqtab(gstab))
  }  
  else
  {
    if(method=="tucker")
    {
      slope1 <- cov(x,xv)/varxv
      slope2 <- cov(y,yv)/varyv
    }
    else if(method=="levine" & internal)
    {
      slope1 <- varx/cov(x,xv)
      slope2 <- vary/cov(y,yv)
    }
    else if(method=="levine" & !internal)
    {
      slope1 <- (varx+cov(x,xv))/(varxv+cov(x,xv))
      slope2 <- (vary+cov(y,yv))/(varyv+cov(y,yv))
    }
    msx <- mx-(yw*slope1*(mxv-myv))
    msy <- my+(w*slope2*(mxv-myv))
    sdsx <- sqrt(varx-
      (yw*(slope1^2)*(varxv-varyv))+
      (w*yw*(slope1^2)*(mxv-myv)^2))
    sdsy <- sqrt(vary+
      (w*(slope2^2)*(varxv-varyv))+
      (w*yw*(slope2^2)*(mxv-myv)^2))
  }  
  out <- list(synthstats=rbind(x=c(mx,sqrt(varx)),y=c(my,sqrt(vary)),
    xv=c(mxv,sqrt(varxv)),yv=c(myv,sqrt(varyv)),xs=c(msx,sdsx),
    ys=c(msy,sdsy)))
  colnames(out$synthstats) <- c("mean","sd")
  if(method=="frequency") 
    out <- c(out,list(synthtab=cbind(scale,xcount=fstab[,2],ycount=gstab[,2])))
  return(out)
}