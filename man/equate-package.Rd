\name{equate-package}
\alias{equate-package}
\docType{package}
\title{
Statistical Methods for Test Score Equating
}
\description{
The equate package contains functions for non-IRT equating under both random groups and nonequivalent groups with anchor test designs. Mean, linear, equipercentile and circle-arc equating are supported, as are methods for univariate and bivariate presmoothing of score distributions. Specific equating methods currently supported include nominal weights, Tucker, Levine observed score, Levine true score, Braun/Holland, frequency estimation, and chained equating.
}
\details{
\tabular{ll}{
Package: \tab equate\cr
Version: \tab 1.1-0\cr
Date: \tab 2010-11-15\cr
Depends: \tab R (>= 2.10.0)\cr
License: \tab GPL (>= 2)\cr
}

Index:
\preformatted{
ACTmath                 ACT Mathematics Test Scores
KBneat                  Test Scores under a NEAT design
descript                Descriptive Statistics for a Frequency Table
equate                  Statistical Equating
freqavg                 Frequency Adjusting and Averaging
freqtab                 Frequency Tabulation
loglinear               Loglinear Smoothing
plot.freqtab            Plotting of Frequency Distributions
px                      Percentile Rank and Cumulative Distribution
                        Functions
se.boot                 Bootstrap Standard Errors of Equating
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{equatevignette} \tab equate vignette (source, pdf)\cr
}
The package vignette provides an introduction to statistical equating and includes descriptions of the supported equating methods and examples. The help page for the main function of the package, \code{\link{equate}}, contains additional examples.
}
\author{
Anthony Albano <tony.d.albano@gmail.com>
}
\keyword{package}
