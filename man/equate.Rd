\name{equate}
\Rdversion{1.1}
\alias{equate}
\title{Linear and Equipercentile Equating}
\description{
This function performs linear and equipercentile equating of test scores under the random groups and nonequivalent groups with anchor test designs (as described in Kolen & Brennan, 2004).
}
\usage{
equate(x, y, scale, type, method = "none", bootse = FALSE,
  internal = TRUE, Ky = max(scale), xv, yv, vscale, w, ...)
}
\arguments{
  \item{x, y}{
vectors of scores for forms X and Y, either as total scores (one per examinee), or total counts (one per score scale point). Total scores are required for the nonequivalent groups design
}
  \item{scale}{
the score scale as a vector
}
  \item{type}{
the type of equating, as the string \code{"mean"}, \code{"linear"} (default), or \code{"equipercentile"}
}
  \item{method}{
string specifying the equating method, where \code{"none"} (default) indicates observed score equating under the random groups design, and \code{"tucker"}, \code{"levine"}, and \code{"frequency"} (i.e., frequency estimation) indicate the respective methods under the nonequivalent groups design
}
  \item{bootse}{
logical indicating whether or not bootstrap standard errors of equating should be returned (default is \code{FALSE})
}
  \item{internal}{
logical indicating whether or not the anchor item scores are included in the total scores (applies only to the Levine method; default is \code{TRUE})
}
  \item{Ky}{
integer, typically indicating the number of items on form Y, or the maximum possible score on Y. This value is used in equipercentile equating for the extremes of the scale (see \code{\link{equate.eq}} for details)
}
  \item{xv, yv}{
vectors of anchor test total scores for form X and Y examinees 
}
  \item{vscale}{
the score scale for the anchor test, as a vector
}
  \item{w}{
value between 0 and 1 specifying the weight applied to form X scores (and implicitly specifying the form Y weight as \code{1 - w})
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
As noted above, for the random groups design \code{x} and \code{y} can be supplied as either total scores, where the length of each will equal the number of examinees, or as counts, where the length of each will equal the length of the score scale with each value indicating the total number of examinees obtaining the corresponding score.

The equating design is implied by the \code{method} argument, where \code{"none"} (default) indicates that no method is needed (because examinees taking forms X and Y are assumed to be randomly equivalent). The Tucker, Levine, and frequency estimation methods are supported for the nonequivalent groups with anchor test design. For these methods only total scores may be used for \code{x} and \code{y} (one per examinee) and  additional arguments are required: \code{xv}, \code{yv}, \code{w} (the weight applied to form X when calculating the synthetic distribution), and \code{vscale}. See \code{\link{synthetic}} for details.

For both mean and linear equating arguments are passed to \code{\link{equate.ln}}, and subsequenty to \code{\link{synthetic}} when the Tucker or Levine method is requested. The Levine method, or more specifically the Levine observed score method, requires the additional argument \code{internal} which is used to estimate the relationship between the anchor and total scores. For equipercentile equating, arguments are passed to \code{\link{equate.eq}}, where smoothing is carried out using the functions \code{\link{freqbump}}, \code{\link{freqavg}} and \code{\link{loglinear}}, each requiring additional smoothing parameters.

Standard errors are returned for the random groups design only. Bootstrap standard errors may be requested, in which case the default is 500 replications, sampling the maximum amount from each score distribution (see \code{\link{se.boot}}).
}
\value{
An object of class \link{equate} is returned, listing the following components (the last 7 are conditional on the equating type, method, and smoothing):
  \item{type}{equating type}
  \item{method}{equating method}
  \item{design}{equating design, inferred from the method}
  \item{stats}{matrix of means and standard deviations for the raw scores}
  \item{freqtab}{frequency table containing the score scale, frequencies, and the cumulative distribution functions for forms X and Y}
  \item{concordance}{conversion table containing scores on X with their form Y equivalents. For the random groups design a column of standard errors is included (see \code{\link{se.ln}} and \code{\link{se.eq}} for details)}
  \item{coefficients}{slope and intercept used in mean and linear equating}
  \item{synthstats}{matrix of means and standard deviations for the anchor and synthetic distributions}
  \item{synthtab}{frequency table for the synthetic distribution, used in frequency estimation equipercentile equating}
  \item{anchortab}{frequency table for the anchor test}
  \item{smoothmethod}{smoothing method}
  \item{smoothout}{a list of smoothing output, the length of which depends on the method chosen (for loglinear smoothing this list will be very long)}
  \item{bootse}{vector of bootstrap standard errors of equating}
}
\references{
Kolen, M. J., & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking}. (2nd ed.), New York: Springer.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{equate}}, \code{\link{equate.ln}}, \code{\link{equate.eq}}, \code{\link{synthetic}}, \code{\link{freqbump}}, \code{\link{freqavg}}, \code{\link{loglinear}}, \code{\link{se.ln}}, \code{\link{se.eq}}, \code{\link{se.boot}}
}
\examples{
#--------------------------------------------------------------
# random groups equating for (1) mean, (2) linear, 
# (3) equipercentile, and (4) equipercentile with 
# loglinear smoothing:

rx <- ACTmath[,2]
ry <- ACTmath[,3]
rscale <- ACTmath[,1]
set.seed(2007)

req1 <- equate(rx,ry,rscale,type="m",bootse=TRUE,reps=100)
req2 <- equate(rx,ry,rscale,type="l",bootse=TRUE,reps=100)
req3 <- equate(rx,ry,rscale,type="e",bootse=TRUE,reps=100)
req4 <- equate(rx,ry,rscale,type="e",bootse=TRUE,reps=100,
  smooth="loglin",degree=3)

# compare equated scores and boostrap standard errors:
cbind(rscale,mean=req1$conc[,2],linear=req2$conc[,2],
  equip=req3$conc[,2],equipS=req4$conc[,2])
cbind(rscale,linear=req2$see,equip=req3$see,equipS=req4$see)

#--------------------------------------------------------------
# nonequivalent groups design for (1) Tucker linear, and
# (2, 3) frequency estimation with weights of 0 and 1

nx <- KBneat$x
ny <- KBneat$y
nscale <- 0:36

neq1 <- equate(nx[,1],ny[,1],nscale,type="Linear",method="Tuck",
  xv=nx[,2],yv=ny[,2],w=1,vscale=0:12)
neq2 <- equate(nx[,1],ny[,1],nscale,type="equip",method="freq",
  xv=nx[,2],yv=ny[,2],w=1,vscale=0:12)
neq3 <- equate(nx[,1],ny[,1],nscale,type="equip",method="freq",
  xv=nx[,2],yv=ny[,2],w=0,vscale=0:12)

# compare equated scores:
cbind(nscale,Tucker=neq1$conc[,2],FEw1=neq2$conc[,2],
  FEw.0=neq3$conc[,2])
}
\keyword{methods}