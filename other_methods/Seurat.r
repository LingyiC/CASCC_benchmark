######### Function
run.Seurat <- function(data.lognormcount, assignfeatures = NULL, resolution = 0.8, PCs= 10, normalize = FALSE, skipPCA = FALSE, nfeatures = 2000){
  res <- list()
  library("Seurat")
  adata <- CreateSeuratObject(counts = data.lognormcount)
  if(is.null(assignfeatures) == TRUE){
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  }else{
    adata@assays$RNA@var.features <- features
  }
  if(normalize == TRUE){
    adata <- NormalizeData(adata, verbose = F)
  }
  
  ## is that ok we skip PCA
  if(skipPCA == FALSE){
    adata <- ScaleData(adata, verbose = F)
    adata <- RunPCA(adata, verbose = F)
    adata <- FindNeighbors(adata, dims = 1:min(ncol(adata@reductions$pca@cell.embeddings), PCs), verbose = F)
    adata <- FindClusters(adata, resolution = resolution, verbose = F)
    adata <- RunUMAP(adata, dim = 1:1:min(ncol(adata@reductions$pca@cell.embeddings), PCs), verbose = F)
  }else {
    adata <- ScaleData(adata, verbose = F)
    # adata <- RunPCA(adata, verbose = F)
    adata <- FindNeighbors(adata, features = adata@assays$RNA@var.features, dim = NULL, verbose = F)
    adata <- FindClusters(adata, resolution = resolution, verbose = F)
    adata <- RunUMAP(adata, features = adata@assays$RNA@var.features, dim = NULL, verbose = F)
  }
  
  clusteringResults <- Idents(adata)
  res$clusteringResults <- clusteringResults
  res$adata <- adata
  return(res)
  
}
#########
