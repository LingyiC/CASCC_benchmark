

##=== RaceID
# Toturial:
# > vignette("RaceID")

# https://github.com/mkrzak/Benchmarking_Clustering_Methods_scRNAseq/blob/master/run_methods/methods/run_RaceID3.R
##===

run.RaceID <- function(data.count, groundTruth, fixedK = FALSE, no_cores = 1, nComp = NULL){
    res <- list()
    library(RaceID)
    # Input to RaceID is a matrix of raw expression values (unique molecular identifiers with or without Poisson correction [@Noise]) with cells as column and genes as rows. This matrix can be provided as a matrix object, a data.frame or a sparse matrix produced by the `Matrix` package.
    sc <- SCseq(data.frame(data.count)) ## To start the analysis, a RaceID single-cell sequencing (SCseq) object is initialized with a count matrix.
    sc <- RaceID::filterdata(sc,mintotal=1)
    # sc <- RaceID::filterdata(sc, minexpr = 5, minnumber = 1) #filters and normalizes data, doesnt work without it
    
    sc <- RaceID::CCcorrect(sc, dimR = TRUE, nComp=nComp) ## Only segerstolpe is FALSE
    sc <- compdist(sc, no_cores = no_cores) # If set to NULL then the number of available cores minus two is used. Default is 1.

    if (fixedK == FALSE) {
        sc = RaceID::clustexp(sc)
    }else if(fixedK == TRUE){
        sc = RaceID::clustexp(sc, sat=FALSE, cln=groundTruth)
    }
    res$sc <- sc
    res$clusteringResults <- as.vector(sc@cluster$kpart)
    detach(package:RaceID)
    return(res)
} 
