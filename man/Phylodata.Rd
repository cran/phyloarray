% Phyloarray data

\name{Phylodata}
\docType{data}
\alias{scans}

\non_function{}

\title{
  Phyloarray example data
}

\format{
  Example data (scans), which is the result from \code{init.data2()}. It
  contains the data of one slide in two colors, washed at nine different
  temperatures.
}

\description{
  }

\author{
  Kurt Sys (\email{kurt.sys@advalvas.be})
}

\seealso{
  \code{\link[phyloarray]{Scandataraw}}

  \code{\link[phyloarray]{init.data}}
}

\keyword{
  datasets
  Phylodata
  scans
}

\examples{
  # load data from \code{scans <- init.data2()}
  data(Phylodata)

  # The attributes of the object \code{scans}
  attributes(scans)

  # The first 50 red values for the different scans
  scans$R[1:10,]

  # The names of the probes
  scans$ID

  # To see the valeus of probe \code{NONEUB_338}
  # All probes are present in triplicate (=> 3 rows)
  scans$R[scans$ID=="NONEUB_338",]
} 


