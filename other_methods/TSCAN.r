##=== TSCAN
# Toturial:
# https://github.com/mkrzak/Benchmarking_Clustering_Methods_scRNAseq/blob/master/run_methods/methods/run_TSCAN.R
##===

run.TSCAN <- function(data.count, groundTruth, fixedK = FALSE, takelog = TRUE){
    res <- list()
    library(TSCAN)

    # mkrzak et al.
    # data <- TSCAN::preprocess(data.lognormcount, takelog = TRUE, logbase = 2, cvcutoff = 0.1) #filter and normalize #default  cvcutoff = 1

    # TUTORIAL https://bioconductor.org/packages/devel/bioc/vignettes/TSCAN/inst/doc/TSCAN.pdf
    # their example data input is raw count
    data <- TSCAN::preprocess(data.count, takelog = takelog, logbase = 2, cvcutoff = 0.1) # default cvcutoff = 1, however, it has error "Error in svd(x, nu = 0, nv = k) : a dimension is zero". 
    # then we follow the mkrzak et al. cvcutoff 

    if(fixedK == TRUE){
      tsc <- TSCAN::exprmclust(data, clusternum = length(levels(groundTruth)), reduce = T)
    }else if(fixedK == FALSE){
      tsc <- TSCAN::exprmclust(data, reduce = T) #PCA #default parameter 
    }
    res$tsc <- tsc
    res$clusteringResults <- as.vector(tsc$clusterid)
    detach(package:TSCAN)
    return(res)
} 
