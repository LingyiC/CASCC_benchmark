.libPaths("/hpc_user_path/rpackages/")


rm(list = ls()); gc()

for (seed.index in 1:5) {
  set.seed(seed.index)
  
  

  library("stringr")
  library(aricode)
  setwd("/hpc_user_path/0215_benchmark_HPC/scanSeeds_")
  source("./RaceID/RaceID.r")

  paths <- list.files("/hpc_user_path/0215_benchmark_HPC/dataset/realGroup_15",
                      pattern=".rds",all.files=FALSE,full.names=T)


  folder_path <- "./RaceID/res"
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
      
    ## RaceID
    start_time <- Sys.time()
    res <- tryCatch({
      run.RaceID(data.count, groundTruth, fixedK = FALSE, no_cores = 1)
    }, error = function(err) {
      run.RaceID(data.count, groundTruth, fixedK = FALSE, no_cores = 1, nComp = 2) ## segerstolpe will do this 
    })
    end_time <- Sys.time()
    
    folder_path <- paste0("./RaceID/res/seed_", seed.index)
    # Check if the folder exists
    if (!file.exists(folder_path)) {
      dir.create(folder_path)
      cat("Folder created:", folder_path, "\n")
    } else {
      cat("Folder already exists:", folder_path, "\n")
    }

    save(res, file = paste0("./RaceID/res/seed_", seed.index, "/RaceID_res_", stringr::str_match(path, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], ".RData"))
    running.time <- c(running.time, as.numeric((end_time-start_time), units = "mins"))

    ## raceID, make 16cell to X16cell, correct it
    vector1 <- res$sc@ndata@Dimnames[[2]]      
    vector2 <- colnames(data.lognormcount)
    for (i in 1:length(vector1)) {
    
      name1 <- vector1[i]
      name2 <- vector2[i]
      
      if(substr(name1, 1, 1) == "X"){
        # Extract the part of the name after the first character
        suffix1 <- substring(name1, 2)
        suffix2 <- substring(name2, 1) 
        # Compare the suffixes
        if (suffix1 == suffix2) {
          vector1[i] <- name2
        }
      }
    }

  
    labels <- res$clusteringResults
    names(labels) <- vector1
    names(groundTruth) <- colnames(data.lognormcount)
    ARI_res <- c(ARI_res, aricode::ARI(as.character(labels),  as.character(groundTruth[names(labels)])))
    AMI_res <- c(AMI_res, aricode::AMI(as.character(labels),  as.character(groundTruth[names(labels)])))
    NMI_res <- c(NMI_res, aricode::NMI(as.character(labels),  as.character(groundTruth[names(labels)])))
  }
  df <- data.frame(dataname = stringr::str_match(paths, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], 
                  ARI_res = ARI_res, 
                  AMI_res = AMI_res, 
                  NMI_res = NMI_res, 
                  running.time = running.time)


  folder_path <- "./RaceID/excel"
  # Check if the folder exists
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
  openxlsx::write.xlsx(df, file = paste0("./RaceID/excel/RaceID_realGroup_15_withoutK_seed_", seed.index, ".xlsx"))


}
