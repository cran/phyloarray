% Print phyloarray-object information

\name{print.phyloarray}

\alias{print.phyloarray}

\title{
  Print phylogeny information
}

\description{
  Prints information of an object of "phyloarray"-class. One should
  however use the generic function \code{print}, which will
  use this function for objects of class "phyloarray".
  
}

\usage{
  print.phyloarray(x, ...)
}

\arguments{
  \item{x}{An object of type \code{phyloarray}}
  \item{...}{Not used}
}

\author{
  Kurt Sys (\email{kurt.sys@advalvas.be})
}

\seealso{
  \code{\link[base]{print}}
}

\keyword{
  }

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  print(scans)
} 
