% Plot meltingcurve

\name{plotmeltingcurve}

\alias{plotmeltingcurve}

\title{
  Plot melting curve
}

\description{
  Make a plot of a melting curve. 
}

\usage{
  plotmeltingcurve(datalist, dye="R", probe=NULL)
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{dye}{The dye of which the background should be plotted, i.e. "R"
	 or "G".}
  \item{probe}{The name of the probe ('Probename') which should be
	 plotted.}
}

\details{
  At present date, only the signal/background ratio is calculated and
  plotted. In the future, other methods may be included.

  The probe-ID's to be plotted are taken from the list of probes. If
  more than one probe has the same 'Probename', all the probes are
  plotted in the figure, using different point character (from type 1 to
  ...). Each next column, i.e. ID of probe of the same 'Probename'
  (match, mismatch, ...) is plotted using a different color, going from
  color 1 to ... See \code{\link[base]{par}} and
  \code{\link[base]{plot}} for more information about characters and
  colors.
  
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

  \code{\link[base]{plot}}
  \code{\link[base]{lines}}
  \code{\link[base]{par}}
}

\keyword{
  melting curve
}

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  # for calculation of histogram cutoffs for good/bad spots
  # scans <- histcutoffs(scans)
  #
  # if background calculation is necessary
  # scans <- calcbackgroud(scans)

  plotmeltingcurve(scans, probe="Bacteria1")
} 
