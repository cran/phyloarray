% Read probes and add to object

\name{getprobes}

\alias{getprobes}

\title{
  Adding probes to object
}

\description{
  Reads the for each "Probename", all match and mismatch probes from a
  \code{.csv}-file. The fields are normally delimited with a comma, and
  the fields surroundend by quotes.
}

\usage{
  getprobes(datalist, file, header=T, sep=",", quote="\"", fill=T, ...)
}

\arguments{
  \item{datalist}{An object of type \code{phyloarray}}
  \item{file}{The file to be read}
  \item{header}{The first row should contain the columns names. The
	 first columns should be named "Probename". The following columns may
	 have any name, and the user should use these names for further
	 analysis.}
  \item{sep}{The delimiter between field. Standard value is a comma.}
  \item{quote}{The quotes around the field values.}
  \item{fill}{If the number of columns is unequal, blank fields are
	 added.}
  \item{...}{Additional parameters for \code{read.csv}}
}

\value{
  An object of class phyloarray is returned, with attribute
  \code{probes} added. This is a dataframe containing for each
  probename, several match and mismatch probe ID's.
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

  \code{\link[base]{read.csv}}
}

\keyword{
  probes
}

\examples{
  # load data \code{Phylodata.RData}, i.e. \code{scans}
  data(Phylodata)

  # read the probes from a \code{.csv}-file
  # scans <- getprobes(scans, file="probes.csv")

  # names of the columns
  names(attr(scans, "probes"))

  # the list of probenames and -ID's
  attributes(scans)$probes
} 
