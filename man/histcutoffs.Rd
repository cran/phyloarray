% Find cutoff values for badspots

\name{histcutoffs}

\alias{histcutoffs}

\title{
  Histogram cutoff values
}

\description{
  Calculate cutoff values for standard deviation/signal intensities for
  badspots. A bad spot has high sd/signal intensity. The cutoff value is
  based on a fivenum/boxplot analysis.
}

\usage{
  histcutoffs(datalist, cutat=1.7)
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray} for which the
	 cutoffs should be calculated.}
  \item{cutat}{A value giving how many times the length of the boxplot
	 the cutoff should be set. This is the \code{coef}-argument of the
	 \code{boxplot.stats} function.}
}

\value{
  The return value is the object datalist of class \code{phyloarray}, but with one
  attribute added, i.e. "cutoff"
}

\details{
  For the function, the standard deviation/signal values are calculated
  and log-transformed. Using the log-transformed values, a boxplot is
  calculated, with the coefficient for outliers at \code{cutat} times
  the length of the box itself.

  Concerning \code{cutat}: If(!) the log-transformed values are normal
  distributed, setting \code{cutat} to 1.0, means having a certainty of
  97.7250\% of no false negatived. 1.35 has 99.3890\% certainty, 1.7 has
  99.8650\% certainty and 2.5 has 99.9968\% certainty of no false
  negatives. It is stressed that these values are only true for normal
  distributions!
}

\note{ 
}

\references{
}

\author{
  Kurt Sys (\email{kurt.sys@advalvas.be})
}

\seealso{
  \code{\link[phyloarray]{Scandataraw}}
  \code{\link[phyloarray]{Phylodata}}

  \code{\link[phyloarray]{init.data}}

  \code{\link[phyloarray]{calcbackground}}

  \code{\link[base]{boxplot.stats}}
  \code{\link[base]{boxplot}}
  \code{\link[base]{fivenum}}
}

\keyword{
  histogram
  boxplot
  cutoff
}

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  # show some histograms
  hist(scans$Rsd[,1]/scans$R[,1], nclass=10000, xlim=c(0,5))
  hist(scans$Gsd[,1]/scans$G[,1], nclass=10000, xlim=c(0,1))

  # calculate cutoff values
  scans <- histcutoffs(scans, cutat=2.5)

  # the cutoff values
  attr(scans, "cutoff")

  # which gives the same as
  attributes(scans)$cutoff
} 
