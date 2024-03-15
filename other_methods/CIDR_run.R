.libPaths("/hpc_user_path/rpackages/")

rm(list = ls()); gc()

for (seed.index in 1:5) {
  set.seed(seed.index)
  
  

  library("stringr")
  library(aricode)
  setwd("/hpc_user_path/0215_benchmark_HPC/scanSeeds")
  source("./CIDR/CIDR.r")

  paths <- list.files("/hpc_user_path/0215_benchmark_HPC/dataset/realGroup_15",
                      pattern=".rds",all.files=FALSE,full.names=T)

  folder_path <- "./CIDR/res"
  # Check if the folder exists
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }

  running.time <- c()
  ARI_res <- c()
  AMI_res <- c()
  NMI_res <- c()
  for(path in paths){
    rds <- readRDS(path)
    data.normcount <- rds@assays[["data"]]@listData[["normcounts"]]
    data.lognormcount <- rds@assays[["data"]]@listData[["logcounts"]]
    groundTruth <- rds@colData@listData[["cell_type1"]]
    groundTruth <- as.factor(groundTruth)

    if(is.null(rds@assays[["data"]]@listData[["counts"]]) == FALSE){
      data.count <- rds@assays[["data"]]@listData[["counts"]]
    }else{
      data.count <- rds@assays[["data"]]@listData[["normcounts"]]
    } 
            
    ## CIDR
    start_time <- Sys.time()
    res <- run.CIDR(data.count, groundTruth, fixedK = FALSE, n_cores = 1)
    end_time <- Sys.time()
    
    folder_path <- paste0("./CIDR/res/seed_", seed.index)
    # Check if the folder exists
    if (!file.exists(folder_path)) {
      dir.create(folder_path)
      cat("Folder created:", folder_path, "\n")
    } else {
      cat("Folder already exists:", folder_path, "\n")
    }

    save(res, file = paste0("./CIDR/res/seed_", seed.index, "/CIDR_res_", stringr::str_match(path, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], ".RData"))
    running.time <- c(running.time, as.numeric((end_time-start_time), units = "mins"))
    ARI_res <- c(ARI_res, aricode::ARI(as.factor(res$clusteringResults), groundTruth))
    AMI_res <- c(AMI_res, aricode::AMI(as.factor(res$clusteringResults), groundTruth))
    NMI_res <- c(NMI_res, aricode::NMI(as.factor(res$clusteringResults), groundTruth))

  }
  df <- data.frame(dataname = stringr::str_match(paths, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], 
                  ARI_res = ARI_res, 
                  AMI_res = AMI_res, 
                  NMI_res = NMI_res, 
                  running.time = running.time)


  folder_path <- "./CIDR/excel"
  # Check if the folder exists
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
  openxlsx::write.xlsx(df, file = paste0("./CIDR/excel/CIDR_realGroup_15_withoutK_seed_", seed.index, ".xlsx"))


}
