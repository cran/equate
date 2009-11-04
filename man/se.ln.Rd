\name{se.ln}
\Rdversion{1.1}
\alias{se.ln}
\title{Standard Error of Linear Equating}
\description{
Computes the standard error of linear equating, under the random groups design, using an equation derived by Braun and Holland (1982).
}
\usage{
se.ln(x, y)
}
\arguments{
  \item{x, y}{
score distributions of class \dQuote{\code{\link{freqtab}}} for forms X and Y
}
}
\details{
This function is called by \code{\link{equate.ln}}, where it is only applicable to the random groups equating design.
}
\value{
Returns a vector of standard errors of equating for each point in the scale
}
\references{
Braun, H. I., & Holland, P. W. (1982). Observed-score test equating: A mathematical analysis of some ETS equating procedures. In P. W. Holland and D. B. Rubin (Eds.), \emph{Test Equating} (pp. 9-49). New York: Academic.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{equate.ln}}
}
\examples{
xscale <- ACTmath[,1]
x <- freqtab(ACTmath[,2],xscale,addclass=TRUE)
y <- freqtab(ACTmath[,3],xscale,addclass=TRUE)
se.ln(x,y)
}
\keyword{univar}