casemod <- function(x,chars=1,upper=TRUE)
{
  if(length(x)>1) stop("'x' must be of length 1")
  n <- nchar(x)
  chars <- chars[1:min(length(chars),n)]
  if(all(chars<=n) & all(chars>0))
  {
    letters <- vector(length=n)
    for(i in 1:n) letters[i] <- substr(x,i,i)
    letters[chars] <- casefold(letters[chars],upper=upper)
    return(paste(letters,collapse=""))
  }
  else(stop("index out of bounds"))
}