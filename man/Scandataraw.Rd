% Raw data-objects

\name{Scandataraw}
\docType{data}

\alias{scan1}
\alias{scan2}
\alias{scan3}

\non_function{}

\title{
  Raw data
}

\description{
  Load the raw example data of one slide (two colors) washed at three
  different temperatures (scan1, scan2, scan3). 
}

\format{3 objects containing raw genepix data}

\author{
  Kurt Sys (\email{kurt.sys@advalvas.be})
}

\seealso{
  \code{\link[phyloarray]{Phylodata}}

  \code{\link[phyloarray]{read.genepix}}
}

\keyword{
  datasets
}

\examples{
  data(Scandataraw)

  # To see the loaded objects
  ls()

  # The header fields of a genepix-file
  scan1$Header

  # The names of the genepix data fields
  names(scan1$Data)
} 
