% Sort data

\name{sortdata}

\alias{sortdata}

\title{
  Sorting data
}

\description{
  Sort data from high to low, based on the data of the first column in
  each data matrix.
}

\usage{
  sortdata(datalist, dye="R", method="subtrbg")
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{dye}{The dye used for sorting data, i.e. "R" or "G".}
  \item{method}{The method for sorting. Only one method is included yet,
  i.e. \code{"subtrbg"}, in which the data is sorted by subtracting the
  background from the signal. More methods will be included in the future.}
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

  # the first 25 rows of sorted scans$R
  sortdata(scans)$R[1:25, ]

  # and their ID's
  sortdata(scans)$ID[1:25]

  # sort scans on green
  scans <- sortdata(scans, dye="G")

  # check probes2probenames
} 
