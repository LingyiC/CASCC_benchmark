remotes::install_version("Seurat", "4.3.0")

rm(list = ls()); gc()

for (seed.index in "default") {
  set.seed(1)


  library("stringr")
  library(aricode)
  library(Seurat)
  setwd("/home/lingyi/Documents/reproduce")


  paths <- list.files("./Datasets/realGroup_15",
                      pattern=".rds",all.files=FALSE,full.names=T)


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
    
    ## CASCC
    start_time <- Sys.time()
    res <- CASCC::run.CASCC(data.lognormcount)
    end_time <- Sys.time()
    save(res, file = paste0("./Evaluation_realgroup15/CASCC_res/oneStepReplicate/CASCC_res_", stringr::str_match(path, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], ".RData"))
    running.time <- c(running.time, as.numeric((end_time-start_time), units = "mins"))
    ARI_res <- c(ARI_res, aricode::ARI(as.factor(res$mainType.output$clusteringResults), groundTruth))
    AMI_res <- c(AMI_res, aricode::AMI(as.factor(res$mainType.output$clusteringResults), groundTruth))
    NMI_res <- c(NMI_res, aricode::NMI(as.factor(res$mainType.output$clusteringResults), groundTruth))

  }
  df <- data.frame(dataname = stringr::str_match(paths, "realGroup_15/\\s*(.*?)\\s*.rds")[, 2], 
                  ARI_res = ARI_res, 
                  AMI_res = AMI_res, 
                  NMI_res = NMI_res, 
                  local.running.time = running.time)
  openxlsx::write.xlsx(df, file = paste0("./Evaluation_realgroup15/CASCC_res/oneStepReplicate/others/oneStep.xlsx"))
}


sink("./Evaluation_realgroup15/CASCC_res/oneStepReplicate/others/session_info.txt")
session_info <- sessionInfo()
print(session_info)
sink()



sink("./Evaluation_realgroup15/CASCC_res/oneStepReplicate/others/system_info.txt")
system_info <- Sys.getenv()
print(system_info)
sink()
