% Read probes and add to object

\name{calcbackground}

\alias{calcbackground}

\title{
  Calculate background
}

\description{
  Background values are calculated based on empty spots. The values for
  non-empty spots are calculated by an interpolation procedure. 
}

\usage{
  calcbackground(datalist, emptyfield="empty")
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{emptyfield}{ID of the spots which are empty, i.e. the spots
	 that are used for background calculation or spots which shouldn't
	 give any signal.}
}

\details{
  The function calculates the background using an kriging algorithm. It
  is important to note that only the spots with \code{ID==emptyfield}
  are used and that badspots, i.e. spots with a sd/signal<cutoff are
  removed from background calculation.

  Other or more flexible interpolation procedures might be included in
  future versions.
}

\value{
  The function returns the same object as given in 'datalist', but the
  \code{Gb} and \code{Rb} values are replaced by the calculated
  background values.
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

  \code{\link[phyloarray]{plotbackground}}

  \code{\link[gstat]{krige}}
}

\keyword{
  background
}

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  # histogram of background values before interpolation procedure
  hist(scans$R[,"40"], nclass=50)

  scans <- calcbackground(scans)

  # histogram of background values after interpolation procedure
  hist(scans$R[,"40"], nclass=50)
} 
