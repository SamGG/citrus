#' Launch web-based interface for configuring and running citrus
#'
#' Launches shiny-based interface to configuring and running Citrus. Creates
#' runCitrus.R that can be used to run analysis in data directory.
#'
#' @param dataDirectory If specified, launches configuration UI with files in
#'   data directory. If \code{NULL}, prompts user to select a single FCS file in
#'   data directory.
#'
#' @param outputDirectory If specified, writes output to specified directory,
#'   otherwise output is written to citrusOutput folder in working directory.
#'
#' @author Robert Bruggner
#' @export
#'
#' @examples
#' # Uncomment to run
#' # citrus.launchUI(file.path(system.file(package = "citrus"),"extdata","example1"))
citrus.launchUI = function(dataDirectory=NULL, outputDirectory=NULL){  
  
  library("shiny")
  library("brew")
  
  if (!is.null(dataDirectory)){
    dataDir <<-dataDirectory
  }
  
  if (is.null(dataDir)) {
    stop("Please specify the directory containing FCS files as an argument or using the global variable dataDir", call. = FALSE)
  }
  
  if (!dir.exists(dataDir)) {
    stop(sprintf("The directory dataDir=\"%s\" does not exist. Please specify the directory containing FCS files as an argument or using the global variable dataDir", dataDir), call. = FALSE)
  }
  
  # dataDir <<- normalizePath(dataDir, .Platform$file.sep)  # absolute path
  
  if (is.null(outputDirectory)) {
    outputDirectory = file.path(dataDir, "citrusOutput")
  }
  if (dir.exists(outputDirectory)) {
    stop(sprintf("Output directory \"%s\" already exists! rename it or erase it first!", outputDirectory), call. = FALSE)
  }
  # output does not exist, thus create it
  dir.create(outputDirectory, showWarnings=F)

  #assign("citrus.outputPath", outputPath, envir = .GlobalEnv)
  citrus.outputPath <<- outputDirectory
  # citrus.outputPath <<- normalizePath(outputDirectory, .Platform$file.sep)  # absolute path
  
  #sapply(list.files(file.path(system.file(package = "citrus"),"shinyGUI","guiFunctions"),pattern=".R",full.names=T),source)
  
  res = tryCatch({
    runApp(appDir=file.path(system.file(package = "citrus"),"shinyGUI"),launch.browser=T)
  }, warning = function(w){
    cat(paste(w,"\n"));
  },error = function(e){
    stop(paste("Unexpected Error:",e))
  }, finally = {
    
  })

  outputPath = citrus.outputPath
  if (runCitrus){
    setwd(outputPath)
    runFile = file.path(outputPath,"runCitrus.R")
    cat(paste("Running Citrus File:",runFile,"\n"))
    logFilePath = file.path(outputPath,"citrusOutput.log")
    cat(paste("Logging output to:",logFilePath,"\n"))
    logFile = .logOn(logFilePath)
    
    source(runFile)
    
    .logOff(logFile=logFile)
  }
  return(paste("Citrus Output in:",outputPath))  
}

#' @export
citrus.getFileParameters = function(fileName,dataDir,...){
  cat(paste0("Reading parameters in ",fileName,"\n"));
  fcsFile = citrus.readFCS(file.path(dataDir,fileName),which.lines=1)
  parameterNames = flowCore::colnames(fcsFile)
  
  # This really should be done as part of citrus.readFCS
  # Same in citrus.readFCSSet
  parameterDescriptions = as.vector(pData(flowCore::parameters(fcsFile))$desc)
  invalidDescriptions = unname(which(sapply(parameterDescriptions,nchar)<3 | is.na(parameterDescriptions)))
  parameterDescriptions[invalidDescriptions] = parameterNames[invalidDescriptions]
  names(parameterNames)=parameterDescriptions
  
  #pnames = as.vector(pData(flowCore::parameters(fcsFile))$desc)
  #pnames[sapply(pnames,nchar)<3] = parameterNames[sapply(pnames,nchar)<3]
  #names(parameterNames)=pnames
  
  return(parameterNames)
}

.logOn = function(filePath,messages=F){
  logFile = file(filePath,open = "wt")
  sink(logFile,split=T)
  if (messages){
    sink(logFile,type="message")
  }
  return(logFile)
}

.logOff = function(logFile,messages=F){
  if (messages){
    sink(type="messsage")
  }
  sink()
  close(logFile)
}
