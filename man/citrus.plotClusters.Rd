% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus.plot.R
\name{citrus.plotClusters}
\alias{citrus.plotClusters}
\title{Plot cluster histograms}
\usage{
citrus.plotClusters(clusterIds, clusterAssignments, citrus.combinedFCSSet,
  clusteringColumns, conditions = NULL, outputFile = NULL, ...)
}
\arguments{
\item{clusterIds}{Vector of cluster IDs to plot}

\item{clusterAssignments}{List containing indicies of cells assigned to each cluster.}

\item{citrus.combinedFCSSet}{Combined FCS data that was clustered.}

\item{clusteringColumns}{Columns for which to plot distributions}

\item{conditions}{Vector of conditions clustering was performed on.}

\item{outputFile}{If not \code{NULL}, plot is written to \code{outputFile}.}

\item{...}{Other parameters (ignored).}
}
\description{
Plot expression of markers in cluster cells relative to all cells
}
\examples{
# Where the data lives
dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")

# Create list of files to be analyzed
fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))

# Read the data 
citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)

# List of columns to be used for clustering
clusteringColumns = c("Red","Blue")

# Cluster data
citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)

# Plot clusters
citrus.plotClusters(clusterIds=c(19998,19997),clusterAssignments=citrus.clustering$clusterMembership,citrus.combinedFCSSet,clusteringColumns)
}
\author{
Robert Bruggner
}
