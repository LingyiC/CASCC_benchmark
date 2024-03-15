
##=== SC3
# Toturial:
# https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html
##===
run.SC3 <- function(data.count, groundTruth, fixedK = FALSE, n_cores = 1){
  res <- list()
  library(SingleCellExperiment)
  library(SC3)
  library(scater)

  groundTruth <- as.factor(groundTruth)
  k.groundTruth <- length(levels(groundTruth))

  ## creat the matrix
  sce <- SingleCellExperiment(
      assays = list(
          counts = as.matrix(data.count),
          ## if the dataset does not provide count matrix
          logcounts = log2(as.matrix(data.count) + 1)    
      ), 
      colData = groundTruth)

  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)


  sce <- SC3::sc3_prepare(sce, gene_filter=FALSE, n_cores = n_cores) 

  if (fixedK == FALSE) {
    sce <- SC3::sc3_estimate_k(sce) 
    cluster.K <- metadata(sce)$sc3$k_estimation
  }else if(fixedK == TRUE){
    cluster.K <- k.groundTruth
  }
  
  sce <- SC3::sc3_calc_dists(sce)
  sce <- SC3::sc3_calc_transfs(sce)
  sce <- SC3::sc3_kmeans(sce, ks = cluster.K)
  sce <- SC3::sc3_calc_consens(sce)  
  
  res$sce <- sce
  res$clusteringResults <- as.factor(colData(sce)[paste0("sc3_",cluster.K , "_clusters")][[1]])
  return(res)
}
