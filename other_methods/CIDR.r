
##=== CIDR
# Tutorial: https://github.com/VCCRI/CIDR
##===
run.CIDR <- function(data.count, groundTruth, n_cores = 1, fixedK = FALSE){
  res <- list()
  library(cidr)
  groundTruth <- as.factor(groundTruth)
  k.groundTruth <- length(levels(groundTruth))

  cid.object <- cidr::scDataConstructor(data.count, tagType = "raw") 
  cid.object <- cidr::determineDropoutCandidates(cid.object)
  cid.object <- cidr::wThreshold(cid.object)
  cid.object <- cidr::scDissim(cid.object, threads = n_cores)
  cid.object <- cidr::scPCA(cid.object, plotPC=FALSE)
  if (fixedK == FALSE) {
            cid.object <- tryCatch({
              cidr::scCluster(object = cid.object, nCluster = NULL, nPC = 10)
            }, error = function(e) {
              message("Error occurred in line1(): ", e$message)
              cidr::scCluster(object = cid.object, nCluster = NULL, nPC = 4)
            
            })   # set nPC = 10 
  }else if(fixedK == TRUE){
    cid.object <- cidr::scCluster(object = cid.object, nCluster = k.groundTruth, nPC = 4) 
  }  
  res$cid.object <- cid.object    
  res$clusteringResults <- cid.object@clusters
  return(res)
}



