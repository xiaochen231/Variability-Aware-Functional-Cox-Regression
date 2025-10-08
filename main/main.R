library(parallel)
files <- list.files(path = "run", full.names = TRUE)

for(files_i in 1:length(files)){
  source(file=files[files_i])
  
  cores <- detectCores(logical = FALSE)
  cl <- makeCluster(cores)
  clusterEvalQ(cl, library(MASS))
  clusterEvalQ(cl, library(grpreg)) 
  clusterEvalQ(cl, library(fda))
  clusterEvalQ(cl, library(survcomp))
  
  
  system.time(
    res1 <- parLapply(cl, 1:100, simulation )
  )
  
  stopCluster(cl)
  
}


