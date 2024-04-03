options(future.globals.maxSize = 500000000 * 1024^2)

library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library("Matrix")
library("Seurat")
library("parallel")
library('randomForest')
library('doParallel')
library('foreach')
library('e1071')
library('SparseM')
## current working directory 
setwd("./")

## Load all functions
source("/Libraries/celltypes_data_gen_t1dvsctrl.R")

# set random seed
set.seed(121) 
#111
##############################################################################################################################
##------------------------------------------------- main resampling methods  ---------------------------------------------##
##############################################################################################################################

freq.matrix = function(train.x, method, ncores, test.x, train.y, test.y)
{
  if(method=="naiveBayes")
  {
    modelnaiveBayes <- naiveBayes(train.x, train.y)
    predictions <- predict(modelnaiveBayes, test.x)
    
    ## confusion matrix
    conf_mat <- confusionMatrix (as.factor(predictions), as.factor(test.y))
    ## accuracy
    acc <- round(conf_mat$overall[[1]]*100,2)
    ## sensitivity
    sens <- round(conf_mat$byClass[[1]]*100,2)
    ## specificity
    spec <- round(conf_mat$byClass[[2]]*100,2)
    ## metrics
    perf <- c(acc, sens, spec)    
    
  }
  return(list(perf))
}

resampling = function(data, y, method, ncores)
{
  ## resampling
  ############## Split data into training and testing ###########################
  s.idx <- sample(1:nrow(data),as.integer(nrow(data)*.7),replace = F)
  
  train.x <- data[s.idx,]
  #train.x <- train.x[1:50,] 
  train.y <- y[s.idx]
  #train.y <- train.y[1:50]
  test.x <-  data[-s.idx,]
  #test.x <- test.x[1:50,]
  test.y <-  y[-s.idx]
  #test.y <- test.y[1:50]
  
  ## xgb
  if(method=="xgb")
  {
    ## Prepare train matrix
    train_matrix <- xgb.DMatrix(data = train.x, label = train.y)
    ## Prepare test matrix
    test_matrix <- xgb.DMatrix(data = test.x, label = test.y)
    
    obj = freq.matrix(train_matrix, test_matrix, method, ncores, test.x, train.y, test.y)
  }
  if(method=="naiveBayes")
  {
    obj = freq.matrix(train.x, method, ncores, test.x, train.y, test.y)
  }  
  # return frequency matrix and coefficients matrix
  return(obj)
}

#################################################################
resampling_results = function(data, y, method, iter, ncores)
{
  # 100 resampling
  # return two lists each of which includes frequency and coefficients
  # for each variables
  #mat0 = list()
  mat1 = matrix(0,nrow=iter,ncol=3)
  #mat2 = list()
  for(i in 1:iter)  ###########################################
  {
    # resampling
    obj0 = resampling(data,y,method,ncores)
    #print(paste0("Iteration ", i, " is completed"))
    #mat0[[i]] = obj0[[1]]
    print(paste0("The output of Iteration ", i, " is ", obj0[[1]]))
    mat1[i,] = obj0[[1]]
    #print(paste0("The cell barcodes are written ", i, " is ", head(obj0[[3]])))
    #mat2[[i]] = obj0[[3]]
  }
  metric.vec = apply(mat1,2,mean)
  sd.vec = apply(mat1,2,sd)
  
  return(list(mat1, metric.vec, sd.vec))
}

main = function(data, y, method, iter, ncores)
{
  results = resampling_results(data, y, method, iter, ncores)
  #genes.vec     = results[[1]]
  mat1 = results[[1]]
  metric.vec = results[[2]]
  sd.vec = results[[3]]
  #cell.df = results[[5]]
  
  return(list(mat1, metric.vec, sd.vec))
}

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
dat <- data.gen()
# dat <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/T1DvsCTL_small_ml.rds")

iter=10
ncores=1

## alpha
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="alpha_naiveBayes_t1dvsctl.rds"
alpha <- dat$alpha
cell.type <- "alpha"
output <- main(alpha[[1]], alpha[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(alpha)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## beta
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="beta_naiveBayes_t1dvsctl.rds"
beta <- dat$beta
cell.type <- "beta"
output <- main(beta[[1]], beta[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(beta)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## acinar
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="acinar_naiveBayes_t1dvsctl.rds"
acinar <- dat$acinar
cell.type <- "acinar"
output <- main(acinar[[1]], acinar[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(acinar)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## ductal
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="ductal_naiveBayes_t1dvsctl.rds"
ductal <- dat$ductal
cell.type <- "ductal"
output <- main(ductal[[1]], ductal[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(ductal)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## delta
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="delta_naiveBayes_t1dvsctl.rds"
delta <- dat$delta
cell.type="delta"
output <- main(delta[[1]], delta[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(delta)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## endothelial
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="endothelial_naiveBayes_t1dvsctl.rds"
endothelial <- dat$endothelial
cell.type <- "endothelial"
output <- main(endothelial[[1]], endothelial[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(endothelial)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## stellates
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="stellates_naiveBayes_t1dvsctl.rds"
stellates <- dat$stellates
cell.type <- "stellates"
output <- main(stellates[[1]], stellates[[2]], "naiveBayes", iter, ncores)
saveRDS(output, filename)
rm(stellates)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## immune
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="immune_naiveBayes_t1dvsctl.rds"
immune <- dat$immune
cell.type = "immune"
output <- main(immune[[1]], immune[[2]], "naiveBayes", iter, ncores) 
saveRDS(output, filename)
rm(immune)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))
