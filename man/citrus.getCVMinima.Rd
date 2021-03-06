% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus.model.R
\name{citrus.getCVMinima}
\alias{citrus.getCVMinima}
\title{Get regularization thresholds of pre-selected cross-validation points}
\usage{
citrus.getCVMinima(modelType, thresholdCVRates, fdrRate = 0.01)
}
\arguments{
\item{modelType}{Method to be used for model-fitting. Valid options are: \code{glmnet},\code{pamr}, and \code{sam}.}

\item{thresholdCVRates}{Matrix of error rates at regularizationThresholds returned by \code{citrus.thresholdCVs.*} function.}

\item{fdrRate}{FDR Maximum used to determine FDR-constrained model regularization threshold.}
}
\value{
List of regularization thresholds and indicies based on pre-selected cross-validation error rates points.
}
\description{
#' Get regularization thresholds of pre-selected cross-validation points and their indicies.
}
\details{
For predictive models (i.e. \code{pamr} or \code{glmnet}), returns indicies of regularization thresholds
producing the minimum cross validation error rate (\code{cv.min}), the simplest model having error within 1
standard error of the minimum (\code{cv.1se}), and the model with the minimum error having an FDR rate < \code{fdrRate} (\code{cv.fdr.constrained})
when possible.
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

# Large enough clusters
largeEnoughClusters = citrus.selectClusters(citrus.clustering)

# Build features
abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)

# List disease group of each sample
labels = factor(rep(c("Healthy","Diseased"),each=10))

# Calculate regularization thresholds
regularizationThresholds = citrus.generateRegularizationThresholds(abundanceFeatures,labels,modelType="pamr",family="classification")

# Calculate CV Error rates
thresholdCVRates = citrus.thresholdCVs.quick("pamr",abundanceFeatures,labels,regularizationThresholds,family="classification") 

# Get pre-selected CV Minima
cvMinima = citrus.getCVMinima("pamr",thresholdCVRates)
}
\author{
Robert Bruggner
}
