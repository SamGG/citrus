import.CITRUS <- function(file,
                          dictionary                = NULL,
                          exclude                   = NULL,
                          bin.width                 = 0.05,
                          minimumClusterSizePercent = 0.05,
                          cluster.selection         = NULL){
  
  message(paste0("Importing ",file))
  if(is.na(file.exists(file)))
    stop(paste0("Error in import.CITRUS: ",file," does not exist"))
  
  env <- new.env()
  load(file,env)
  citrus.combinedFCSSet <- env[["citrus.combinedFCSSet"]]
  citrus.foldClustering <- env[["citrus.foldClustering"]]
  citrus.clustering     <- citrus.foldClustering$allClustering
  
  # from citrus.selectClusters.minimumClusterSize
  browser()
  citrus.allCluster   <- citrus.foldClustering$allClustering
  cluster.sizes       <- sapply(citrus.clustering$clusterMembership,length)
  size.min            <- (length(cluster.sizes)+1) * minimumClusterSizePercent
  largeEnoughClusters <- which(cluster.sizes >= size.min)
  # from citrus.selectClusters.minimumClusterSize
  
  if(!is.null(cluster.selection)){
    largeEnoughClusters <- largeEnoughClusters[largeEnoughClusters%in%cluster.selection]
    if(length(largeEnoughClusters) == 0)
      stop("Selected clusters do not have enough cell")
  }
  data <- c()
  for(clusterId in largeEnoughClusters){
    
    # from citrus.exportCluster - begin
    clusterData = citrus.combinedFCSSet$data[citrus.clustering$clusterMembership[[clusterId]],]
    if (!is.null(citrus.combinedFCSSet$scaleColumns)) {
      clusterData[, citrus.combinedFCSSet$scaleColumns] = t((t(clusterData[,
                                                                           citrus.combinedFCSSet$scaleColumns]) * citrus.combinedFCSSet$scaleColumns.SD) +
                                                              citrus.combinedFCSSet$scaleColumns.mean)
    }
    if (!is.null(citrus.combinedFCSSet$transformColumns)) {
      clusterData[, citrus.combinedFCSSet$transformColumns] = sinh(clusterData[,
                                                                               citrus.combinedFCSSet$transformColumns]) * citrus.combinedFCSSet$transformCofactor
    }
    # from citrus.exportCluster - end
    
    data <- rbind(data,clusterData)
  }
  ##data <- citrus.combinedFCSSet$data
  
  
  if(!is.null(dictionary)){
    colnames(data) <- rename.markers(colnames(data),dictionary)
  }
  if(!is.null(exclude)){
    data <- exclude.markers(data,exclude)
  }
  markers        <- colnames(data)
  
  cluster <- as.CLUSTER(data,cluster = "fileEventNumber",bin.width = bin.width)
  cluster@name <- basename(file)
  
  
  return(cluster)
}

