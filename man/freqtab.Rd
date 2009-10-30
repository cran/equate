\name{freqtab}
\Rdversion{1.1}
\alias{freqtab}
\title{Frequency Tabulation}
\description{
This function creates a frequency table for a vector of scores.
}
\usage{
freqtab(x, scale)
}
\arguments{
  \item{x}{
vector of scores either of length greater than \code{length(scale)}, where each value represents an examinee total score, or of length equal to \code{length(scale)}, where each value is already a score frequency (see below for details)
}
  \item{scale}{
the score scale, as a vector
}
}
\details{
If \code{x} is a vector of total scores, values are factored with \code{levels = scale}. This ensures that unobtained scores (i.e., scores not in \code{x}) will still be represented in the tabulation. If \code{x} is a vector of frequencies, with one value for each scale point provided in \code{scale}, the function serves to simply change the class of \code{x} to \code{\link{freqtab}}.
}
\value{
Returns a 2-column matrix, where column 1 contains the score scale and column 2 the frequency for each scale value
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{table}}
}
\examples{
set.seed(2005)
x1 <- round(rnorm(100,100,10))
x1scale <- 1:100
freqtab(x1,x1scale)

x2 <- round(rnorm(1000,100,10))
x2scale <- 70:130
freqtab(x2,x2scale)

x3scale <- ACTmath[,1]
x3 <- ACTmath[,2]
freqtab(x3,x3scale)
}
\keyword{distribution}
