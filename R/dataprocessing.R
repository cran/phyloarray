require(gstat) || stop("Can't find package gstat")



print.phyloarray <- function(x, ...) {
  cat(paste("    n: ", length(x$ID), "\n", sep=""))
  cat("    columns: "); cat(colnames(x$R)); cat("\n\n")
  cat(paste("  Class: ", attr(x, "class"), "\n\n", sep=""))
  cat(paste("  Probes:\n"))
  cat(paste("    n:", length(attr(x, "probes")[,1])), "\n", sep="")
  cat("    columns: "); cat(names(attr(x, "probes"))); cat("\n\n")
}



probes2probenames <- function(datalist, probe) {
  phyl <- attr(datalist, "probes")
  vals <- NULL
  for (i in probe)
    for (p in 2:length(names(phyl)))
      vals <- c(vals, levels(phyl$Probename)[phyl$Probename[phyl[[p]]==i]])
  return(levels(factor(vals)))
}





sortdata <- function(datalist, dye="R", method="subtrbg") {
  vals <- NULL
  tmp.datalist <- datalist[[dye]] - datalist[[paste(dye, "b", sep="")]]
  for (p in 1:dim(tmp.datalist)[1])
    vals <- c(vals, max(tmp.datalist[p,]))
  for (f in names(datalist))
    if (length(dim(datalist[[f]]))>1)
      datalist[[f]] <- datalist[[f]][sort.list(vals, decreasing=TRUE),]
    else
      datalist[[f]] <- datalist[[f]][sort.list(vals, decreasing=TRUE)]
  return(datalist)
}





plotmeltingcurve <- function(datalist, dye="R", probe=NULL)
{
  phyl <- attr(datalist, "probes")
  if (length(probe)==0)
    probe <- phyl$Probename[1]
  vals <- list()
  maximum <- 0
  for (m in names(phyl)[2:dim(phyl)[2]]) {
    vals[[m]] <- list()
    tmp.vals <- NULL
    for (p in levels(phyl[[m]])[phyl[[m]][phyl$Probename==probe]]) {
      tmp.vals <- datalist[[dye]][datalist$ID==p,]/datalist[[paste(dye,"b",sep="")]][datalist$ID==p,]
      maximum <- max(maximum, tmp.vals)
      if (dim(tmp.vals)[1] > 1)
        vals[[m]][[p]] <- data.frame(rbind(apply(tmp.vals, 2, mean), apply(tmp.vals,2,sd)), row.names=c("Mean", "SD"))
      else
        vals[[m]][[p]] <- tmp.vals
      colnames(vals[[m]][[p]]) <- colnames(tmp.vals)
    }
  }
  plot(as.numeric(vals[[1]][[1]]["Mean",]) ~ as.numeric(colnames(vals[[1]][[1]])), ylim=c(0, maximum), xlab="Temperature", ylab="Signal/Background")
  if (dim(vals[[1]][[1]])[1] > 1)
    arrows(as.numeric(colnames(vals[[1]][[1]])),
           as.numeric(vals[[1]][[1]]["Mean",]) - as.numeric(vals[[1]][[1]]["SD",]),
           as.numeric(colnames(vals[[1]][[1]])),
           as.numeric(vals[[1]][[1]]["Mean",]) + as.numeric(vals[[1]][[1]]["SD",]),
           code=3,
           angle=90,
           length=0.05)
  character <- 1
  color <- 1
  for (m in names(vals)) {
    for (p in names(vals[[m]])) {
      if (length(vals[[m]][[p]]) > 0) {
        lines(as.numeric(vals[[m]][[p]]["Mean",]) ~ as.numeric(colnames(vals[[m]][[p]])), pch=character, col=color, type="p")
        if (dim(vals[[m]][[p]])[1] > 1)
          arrows(as.numeric(colnames(vals[[1]][[1]])),
                 as.numeric(vals[[m]][[p]]["Mean",]) - as.numeric(vals[[m]][[p]]["SD",]),
                 as.numeric(colnames(vals[[m]][[p]])),
                 as.numeric(vals[[m]][[p]]["Mean",]) + as.numeric(vals[[m]][[p]]["SD",]),
                 code=3,
                 angle=90,
                 length=0.05,
                 col=color)
        character <- character+1
      }
    }
    character <- 1
    color <- color+1
  }
}





plotbackground <- function(datalist, dye="R", scan=1)
{
  filled.contour(x=1:(max(as.numeric(datalist$X)) - min(as.numeric(datalist$X)) + 1),
                 y=1:(max(as.numeric(datalist$Y)) - min(as.numeric(datalist$Y)) + 1),
                 z=matrix(datalist[[paste(dye, "b", sep="")]][,scan], ncol=(max(as.numeric(datalist$Y)) - min(as.numeric(datalist$Y)) + 1)),
                 nlevels=100)
}






calcbackground <- function(datalist, emptyfield="empty")
{
  for (dye in c("R", "G")) {
    for (slide in 1:length(colnames(datalist[[dye]]))) {
      emptyspot <- (datalist$ID==emptyfield) & (datalist[[paste(dye, "sd", sep="")]][,slide]/datalist[[dye]][,slide] < attr(datalist, "cutoff")[dye, slide])
      v <- list()
      v$x <- as.numeric(datalist$X[emptyspot])
      v$y <- as.numeric(datalist$Y[emptyspot])
      v$z <- as.numeric(datalist[[dye]][,slide][emptyspot])
      grid <- data.frame(x=as.numeric(datalist$X), y=as.numeric(datalist$Y))
      tmp <- krige(z ~ 1, ~x+y, newd=grid, data=v, block=c(1,1))
      datalist[[paste(dye, "b", sep="")]][,slide] <- c(tmp$var1.pred)
    }
  }
  return(datalist)
}





histcutoffs <- function(datalist,
                        cutat=1.7)
{
  d.cutoff <- NULL
  for (dye in c("R", "G")) {
    tmp.cutoff <- NULL
    tmp.datalist <- log(datalist[[paste(dye, "sd", sep="")]]/datalist[[dye]])
    for (slide in 1:length(colnames(datalist[[dye]]))) {
      tmp.data <- tmp.datalist[,slide]
      tmp.boxstats <- boxplot.stats(tmp.data, coef=cutat)$stats
                                        # range=1: +1 time length of box (relative to box edges)
                                        # logic: if normal: 1 sd, 84% of the points
                                        #                   0.6745 sd, 75% interval (or 25% higher than mean)
                                        # for:   2 times sd, 97.7250% certainty no false 'negatives' -> pnorm(2)
                                        #        2/qnorm(0.75) is about +1.0 times box length
                                        # for:   2.5 times sd, 99.3890% certainty no false 'negatives'
                                        #        2.5/qnorm(0.75) is about +1.35 times box length
                                        # for:   3 times sd, 99.8650% certainty no false 'negatives'
                                        #        3/qnorm(0.75) is about +1.7 times box length
                                        # for:   4 times sd, 99.9968% certainty no false 'negatives'
                                        #        4/qnorm(0.75) is about +2.5 times box length
                                        # distribution might not be normal, but still...
      tmp.cutoff <- c(tmp.cutoff,
                      exp(tmp.boxstats[5]))
    }
    d.cutoff <- rbind(d.cutoff, tmp.cutoff)
  }
  rownames(d.cutoff) <- c("R", "G")
  attr(datalist,"cutoff") <- d.cutoff
  return(datalist)
}





getprobes <- function(datalist, file, header=T, sep=",", quote="\"", fill=T,...)
{
  tmp.frame <- read.csv(file, header=header, sep=sep, quote=quote, fill=fill, ...)
  attr(datalist, "probes") <- tmp.frame
  return(datalist)
}






init.data2 <- function ()
{
  cat("Are you creating a new data matrix or adding new array data\n")
  cat("to a prexisting data matrix? \n")
  cat("Enter \"n\" for creating  and \"a\" for adding new array data: ")
  new.n <- readline()
  if (new.n == "a") {
    cat("Enter the name of the existing data matrix: ")
    oname <- readline()
  }
  cat("Do the names of all your datasets have the following format: \n")
  cat("prefix1, prefix2, prefix3?, ... Here prefix can be any name, \n")
  cat("but the suffixes must be integers 1,2, ..., # of arrays. \n")
  cat("Enter \"y\" for yes, \"n\" for no: ")
  b.n <- readline()
  if (b.n == "y") {
    cat("Enter the prefix:")
    prefixname <- readline()
    cat("Enter the number of arrays to be processed:")
    n <- readline()
    n <- as.integer(n)
    dname <- paste(prefixname, 1:n, sep = "")
  }
  else if (b.n == "n") {
    cat("Enter the number of arrays to be processed:")
    n <- as.integer(readline())
    dname <- rep(0, n)
    for (i in 1:n) {
      cat(paste("Enter the name of your ", i, "th dataset:"))
      dname[i] <- readline()
    }
  }
  cat("Enter the name of Cy3 raw data: ")
  name.G <- readline()
  cat("Enter the name of Cy3 background: ")
  name.Gb <- readline()
  cat("Enter the name of Cy3 intensity standard deviation: ")
  name.Gsd <- readline()
  cat("Enter the name of Cy5 raw data: ")
  name.R <- readline()
  cat("Enter the name of Cy5 background: ")
  name.Rb <- readline()
  cat("Enter the name of Cy5 intensity standard deviation: ")
  name.Rsd <- readline()
  if (new.n == "a") {
    res <- eval(as.name(oname))
    action <- "updating"
  }
  else {
    res <- list(R = NULL, G = NULL, Rb = NULL, Gb = NULL, Rsd=NULL, Gsd=NULL, ID=NULL)
    cat("Enter the name of the ID/probenames field: ")
    name.ID <- readline()
    res$ID <- as.character(eval(as.name(dname[1]))$Data[, name.ID])
    cat("Enter the X values for probe positions: ")
    name.X <- readline()
    res$X <- as.character(eval(as.name(dname[1]))$Data[, name.X])
    cat("Enter the Y values for probe positions: ")
    name.Y <- readline()
    res$Y <- as.character(eval(as.name(dname[1]))$Data[, name.Y])
    action <- "creating"
  }
  for (i in c("R", "G", "Rb", "Gb", "Rsd", "Gsd"))
    vstr.colnames <- colnames(res[[i]])
  for (i in 1:n) {
    vstr.colnames <- c(vstr.colnames, eval(as.name(dname[i]))$Header$Temperature)
    tmp <- eval(as.name(dname[i]))$Data[, c(name.R, name.G, name.Rb,
                                            name.Gb, name.Rsd, name.Gsd)]
    res$R <- cbind(res$R, as.numeric(as.vector(tmp[, 1])))
    res$G <- cbind(res$G, as.numeric(as.vector(tmp[, 2])))
    res$Rb <- cbind(res$Rb, as.numeric(as.vector(tmp[, 3])))
    res$Gb <- cbind(res$Gb, as.numeric(as.vector(tmp[, 4])))
    res$Rsd <- cbind(res$Rsd, as.numeric(as.vector(tmp[, 5])))
    res$Gsd <- cbind(res$Gsd, as.numeric(as.vector(tmp[, 6])))
  }
  for (i in c("R", "G", "Rb", "Gb", "Rsd", "Gsd"))
    colnames(res[[i]]) <- vstr.colnames
  cat(paste("Finished", action, "the dataset.\n", sep = " "))
  attr(res, "class") <- "phyloarray"
  return(res)
}





read.genepix2 <- function(file,
                          temperature=42,
                          buffer=data.frame(na=0.97, f=35))
{
  lstr.lines <- readLines(file)
  data <- list()
  data$Header <- list()
  i.line <- 0
                                        # Read header first...
  repeat  {
    i.line <- i.line + 1
    if(substr(lstr.lines[i.line],1,1) == "\"" && length(lstr.lines[i.line]) > 0)
        lstr.lines[i.line] <- eval(parse(text=lstr.lines[i.line]))
    if((regexpr("Block", lstr.lines[i.line+1]) > 0) & (regexpr("Name", lstr.lines[i.line+1]) > 0))
        break;
    tmp.str <- strsplit(lstr.lines[i.line], "=")[[1]]
    if (length(tmp.str) > 1)
      data$Header[[ tmp.str[1] ]] <- strsplit(tmp.str[2:length(tmp.str)], "\t")[[1]]
    else
      data$Header[[ i.line ]] <- strsplit(tmp.str, "\t")[[1]]
  }

  data$Header$Temperature <- temperature
  data$Header$Buffer <- buffer
  i.line <- i.line
                                        # ... and now, read data
  data$Data <- read.table(file, skip=i.line, header=T)
  gpname <- colnames(data$Data)
  gpname <- gsub("F635.", "R", gpname, ignore.case=F)
  gpname <- gsub("F532.", "G", gpname, ignore.case=F)
  gpname <- gsub("B635.", "Rb", gpname, ignore.case=F)
  gpname <- gsub("B532.", "Gb", gpname, ignore.case=F)
  gpname <- gsub("Median", "med", gpname, ignore.case=F)
  gpname <- gsub("Mean", "mean", gpname, ignore.case=F)
  gpname <- gsub("SD", "SD", gpname, ignore.case=F)
  colnames(data$Data) <- gpname
  return(data)
}
