% Find probenames

\name{probes2probenames}

\alias{probes2probenames}

\title{
  Find probenames
}

\description{
  Search probenames if a probe-ID is known.
}

\usage{
  probes2probenames(datalist, probe)
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{probe}{ID of the probe of which the probename is searched.}
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

  # sort scans on green
  scans <- sortdata(scans, dye="G")

  plotprobes <- probes2probenames(scans, scans$ID[1:25])

  for (p in plotprobes)
    plotmeltingcurve(scans, dye="G", p)
} 
