\name{casemod}
\Rdversion{1.1}
\alias{casemod}
\title{Lettercase Modification}
\description{
This function changes the case of specified letters within a string
}
\usage{
casemod(x, chars = 1, upper = TRUE)
}
\arguments{
  \item{x}{
string
}
  \item{chars}{
vector of integers indexing the characters to be modified (see below for details)
}
  \item{upper}{
logical. If \code{TRUE} (default), specified characters will be set to uppercase, otherwise lowercase
}
}
\details{
By default the first letter is capitalized. \code{chars} is truncated to be of length \code{length(x)} and may only contain values between 1 and \code{nchar(x)}
}
\value{
Returns a string
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{casefold}}, \code{\link{substr}}
}
\examples{
casemod("latex")
casemod("latex",c(1,3,5))
}
\keyword{character}