% Read genepix-files

\name{read.genepix}

\alias{read.genepix2}

\title{
  Reading genepix-file
}

\description{
  Read a genepix-file ('.gpr') and stores result in a list containing
  'header' and 'data'. 
}

\usage{
  read.genepix2(file, temperature=42, buffer=data.frame(na=0.97, f=35))
}

\arguments{
  \item{file}{The name of the file to be read. The extension should be
	 given.}
  \item{temperature}{The temperature of hybridization (or washing)
	 step. Since there's a focus on melting curves, this is important and
	 the value given here will be added to the header.}
  \item{buffer}{Buffer contents used during hybridization (or
	 washing). This is not used anywhere yet, but it might e.g. to
	 calculate the theoretical melting/annealing temperature based on the
	 sequence and the buffer content. 'na' is the salt (Na)--content and
	 'f' the formamide concentration (in percent).}
}

\value{
  The function returns a list containing:
  \item{Header}{A list with the header fields of the \code{.gpr}-file, with
	 additional information concerning the hybridization temperature and
	 buffer contents (given as arguments to the function.}
  \item{Data}{A data matrix with the data in a \code{.gpr}-file. This
	 includes median and mean values for fore- and background, as well as
	 probe ID's, standard deviations, X and Y	values etc.}
}

\note{
  Since there is a function read.genepix in \code{sma}, this function has
  been called \code{read.genepix2}. It must be noted that for new versions of
  genepix, the function in the \code{sma}-package doesn't work (since more
  columns are added). This function works for both, but I have really no
  idea how the future will be... Anyway I tried to make it as general as
  possible.

  Data from other image processing software (e.g. Spot) may be imported
  using the package \code{sma}. In that case, one should know the names
  of the fields of the red/green signals etc.
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
}

\keyword{
  genepix file
  raw data
}

\examples{
  # There are \code{.gpr}-files in the data-section of this package
  #
  # scan1 <- read.genepix2("<file_scan.at.22degrees>", temperature=22)
  # ...
  # scan3 <- read.genepix2("<file_scan.at.58degrees>", temperature=58)
  #
  # 
  # or, e.g.:
  #
  # scan <- 1
  # for (t in c(22+0:2*18)) {
  #   assign(paste("scan", slide, sep=""),
  #          read.genepix2(paste("wash", t, ".gpr", sep=""), temperature=t))
  #   scan <- scan+1
  # }
}
