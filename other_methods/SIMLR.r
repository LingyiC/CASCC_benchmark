
##=== SIMLR 
# Tutorial: # https://github.com/BatzoglouLabSU/SIMLR/blob/SIMLR/R_main_demo_SIMLR.R
# The default value of cores.ratio in SIMLR is 1, which means that all available CPU cores will be used for parallel processing by default. You can adjust this value based on the specific resources available on your system and the demands of your analysis.
##===
run.SIMLR <- function(data.lognormcount, groundTruth, fixedK = FALSE, cores.ratio = 1){
  res <- list()
  

  if (fixedK == FALSE) {
    k.predict <- SIMLR::SIMLR_Estimate_Number_of_Clusters(data.lognormcount, NUMC = 2:15)
    cluster.k <- seq(2,15)[which.min(k.predict$K1)]
    res$k.predict <- k.predict
    res$cluster.k <- cluster.k
  }else if(fixedK == TRUE){
    k.groundTruth <- length(levels(as.factor(groundTruth)))
    cluster.k <- k.groundTruth

  }

  res$simlr.obj <- SIMLR::SIMLR(X=data.lognormcount, c=cluster.k, normalize = FALSE, cores.ratio = cores.ratio) # `normalize` should I normalize the input data (default = F)? 
  res$clusteringResults <- as.factor(res$simlr.obj$y$cluster)
  return(res)
}
