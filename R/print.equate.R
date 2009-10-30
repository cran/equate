print.equate <- function(x,...)
{
  design <- casemod(x$design)
  if(x$method=="none") method <- casemod(x$type)
  else method <- paste(casemod(x$method),x$type)
  
  cat("\n",method," equating: ",design,"\n",sep="")
  cat("\nSummary Statistics:\n")
  if(design=="random groups") print.default(x$stats)
  else print.default(rbind(x$stats,x$synthstats),quote=FALSE)
  if(method!="Equipercentile")
  {
    cat("\nCoefficients:\n")
    print.default(x$coef,quote=FALSE)
  }
  cat("\nConcordance Table:\n")
  print.default(x$concord,quote=FALSE)
  invisible(x)
}