% Read important fields from raw data-objects

\name{init.data}

\alias{init.data2}

\title{
  Initialize data
}

\description{
  Read the important fields of the several raw data-objects and store
  them in a new object (with class "phyloarray"). This function also
  reads the temperatures from the headers.
}

\usage{
  init.data2()
}

\arguments{
}

\value{
  The function returns an object of class \code{phyloarray}. It is a list
  with values R (red signal), G (green signal), Rb (red background
  signal), Gb (green background signal), Rsd (standard deviation of red
  signal) and Gsd (standard deviation of green signal). The standard
  deviation values are used to mark bad spots (based on sd/signal). The
  list contains also X- and Y-values of the spots (to calculate
  background, if necessary) and the ID's of the probes.
  Each field R, G, Rb, Gb, Rsd and Gsd is a datamatrix with in the rows
  the values for the different probes and in the columns the different
  temperatures. To a larger extent, the columns may be other factors to
  (such as time-related values), but they are designed to be used for
  making melting curves and analysis on the melting curves.
}

\note{
  Since there is a function \code{init.data} in \code{sma}, this function has
  been calles \code{init.data2}. It is largely the same function, but some
  fields are added to make the phylogenetic analysis of melting curves
  more easy.

  Some examples of basic handling of an object of class phyloarray is in
  the help-file \code{Phylodata}
}

\references{
}

\author{
  Kurt Sys (\email{kurt.sys@advalvas.be})
}

\seealso{
  \code{\link[phyloarray]{Scandataraw}}
  \code{\link[phyloarray]{Phylodata}}

  \code{\link[phyloarray]{read.genepix}}
  
  \code{\link[phyloarray]{getprobes}}
}

\keyword{
  initialize data
  phyloarray
}

\examples{
  # load the saved raw data from this package:
  # data(Scandataraw)

  # scans <- init.data2()

  # Are you creating a new data matrix or adding new array data
  # to a prexisting data matrix? 
  # Enter "n" for creating  and "a" for adding new array data: n
  # Do the names of all your datasets have the following format: 
  # prefix1, prefix2, prefix3?, ... Here prefix can be any name, 
  # but the suffixes must be integers 1,2, ..., # of arrays. 
  # Enter "y" for yes, "n" for no: y
  # Enter the prefix:scan
  # Enter the number of arrays to be processed:3
  # Enter the name of Cy3 raw data: Gmed
  # Enter the name of Cy3 background: Gbmed
  # Enter the name of Cy3 intensity standard deviation: GSD
  # Enter the name of Cy5 raw data: Rmed
  # Enter the name of Cy5 background: Rbmed
  # Enter the name of Cy5 intensity standard deviation: RSD
  # Enter the name of the ID/probenames field: ID
  # Enter the X values for probe positions: Column
  # Enter the Y values for probe positions: Row
  # Finished creating the dataset.



  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  names(scans)

  attributes(scans)

} 


