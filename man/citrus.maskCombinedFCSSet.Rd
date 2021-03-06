% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus.util.R
\name{citrus.maskCombinedFCSSet}
\alias{citrus.maskCombinedFCSSet}
\title{Masks a citrus.combinedFCSSet}
\usage{
citrus.maskCombinedFCSSet(citrus.combinedFCSSet, fileIds)
}
\arguments{
\item{citrus.combinedFCSSet}{A citrus.combinedFCSSet object.}

\item{fileIds}{Vector of file IDs for which to retain data.}
}
\value{
A \code{citrus.combinedFCSSet} object.
}
\description{
Masks a citrus.combinedFCSset to include data from a subset of file ids.
}
\examples{
# Where the data lives
dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")

# Create list of files to be analyzed
fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))

# Read citrus.combinedFCSSet
citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)

# Mask
maskedFCSSet = citrus.maskCombinedFCSSet(citrus.combinedFCSSet,fileIds=1:3)

# Check - should be 1,2,3
unique(maskedFCSSet$data[,"fileId"])
}
\seealso{
\code{\link{citrus.readFCSSet}}
}
\author{
Robert Bruggner
}
