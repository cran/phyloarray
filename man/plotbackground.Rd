% Plot background

\name{plotbackground}

\alias{plotbackground}

\title{
  Plot background
}

\description{
  Make a plot of the background of a certain color and scan.
}

\usage{
  plotbackground(datalist, dye="R", scan=1)
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{dye}{The dye of which the background should be plotted, i.e. "R"
	 or "G".}
  \item{scan}{The scan of which the background should be plotted. This
	 may be a number (the number of the column) or a string (name of the column).}
}

\value{
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

  \code{\link[phyloarray]{calcbackground}}

  \code{\link[base]{filled.contour}}
}

\keyword{
  background calculation
}

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  # plot background for red, temperature=58°C
  plotbackground(scans, dye="R", scan="58")

  # this is the same
  plotbackground(scans, scan="58")

  # and this too
  plotbackground(scans, scan=3)
} 
