library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library("Matrix")
library("Seurat")
library("parallel")
## current working directory 
setwd("./")

## Load all funcitons
source("/Libraries/celltypes_data_gen_third.R")

# set random seed
set.seed(999)
##############################################################################################################################
##------------------------------------------------- main resampling methods  -----------------------------------------------##
##############################################################################################################################

freq.matrix = function(train_matrix, test_matrix, method, ncores, train.y, test.y, test.x)
{
  ## XGB
  if(method=="xgb")
  {
    best_param = list()
    best_seednumber = 9999
    best_logloss = Inf
    best_logloss_index = 0
    
    for (iter in 1:10) {
      param <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    max_depth = sample(3:7, 1),
                    eta = runif(1, .01, .3),
                    gamma = runif(1, 0.0, 0.2), 
                    subsample = runif(1, .5, .8),
                    colsample_bytree = runif(1, .5, .8), 
                    min_child_weight = sample(1:30, 1)
                    # max_delta_step = sample(1:10, 1)
      )
      cv.nround = 200
      cv.nfold = 5
      seed.number = sample.int(10000, 1)[[1]]
      set.seed(seed.number)
      mdcv <- xgb.cv(data=train_matrix, params = param, nthread=ncores, 
                     nfold=cv.nfold, nrounds=cv.nround,
                     verbose = T, print_every_n = 50, 
                     early_stopping_rounds=100, maximize=FALSE)
      min_logloss = min(mdcv$evaluation_log[, test_logloss_mean])
      min_logloss_index = which.min(mdcv$evaluation_log[, test_logloss_mean])
      # min_logloss_index <- mdcv$best_iteration
      if (min_logloss < best_logloss) {
        best_logloss = min_logloss
        best_logloss_index = min_logloss_index
        best_seednumber = seed.number
        best_param = param
      }
    }
    
    nround = best_logloss_index
    set.seed(best_seednumber)
    bst_model <- xgb.train(data=train_matrix, params=best_param, nrounds=nround, nthread=ncores)
    
    # Predict hold-out test set
    test_pred <- predict(bst_model, newdata = test_matrix)
    test_pred <- ifelse (test_pred > 0.5,1,0)
    names(test_pred) <- rownames(test.x)
  }
  return(test_pred)
}

resampling = function(dat, method, ncores)
{
  ## resampling
  ##### ----------------------Split data------------------------------#####
  ## Train data
  dat_T1DvsCtrl <- subset(dat, subset = disease_state != "AAB")
  train.x <-GetAssayData(dat_T1DvsCtrl, assay="SCT",slot="data")
  train.x <- t(train.x)
  train.y <- dat_T1DvsCtrl$class_label_num
  train.y <- ifelse(train.y == 2, 1, train.y)
  
  ## Test data 
  dat_AAB <- subset(dat, subset = disease_state == "AAB")
  test.x <-GetAssayData(dat_AAB, assay="SCT",slot="data")
  test.x <- t(test.x)
  test.y <- dat_AAB$class_label_num
  test.y <- ifelse(test.y == 1, 0, test.y)
  
  ## Prepare train matrix
  train_matrix <- xgb.DMatrix(data = train.x, label = train.y)
  ## Prepare test matrix
  test_matrix <- xgb.DMatrix(data = test.x, label = test.y)
  
  ## xgb
  if(method=="xgb")
  {
    obj = freq.matrix(train_matrix, test_matrix, method, ncores, train.y, test.y, test.x)
  }
  # return frequency matrix and coefficients matrix
  return(obj)
}

##################################################################
## Main function
main = function(dat, method, ncores)
{
  results = resampling(dat, method,ncores)
  return(results)
}

dat <- data.gen()

ncores=30
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="alpha_xgb_third.rds"
alpha <- dat$alpha
output <- main(alpha, "xgb", ncores) 
saveRDS(output, filename)
rm(alpha)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="beta_xgb_third.rds"
beta <- dat$beta
output <- main(beta, "xgb", ncores) 
saveRDS(output, filename)
rm(beta)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="acinar_xgb_third.rds"
acinar <- dat$acinar
output <- main(acinar, "xgb", ncores) 
saveRDS(output, filename)
rm(acinar)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="ductal_xgb_third.rds"
ductal <- dat$ductal
output <- main(ductal, "xgb", ncores) 
saveRDS(output, filename)
rm(ductal)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="immune_xgb_third.rds"
immune <- dat$immune
output <- main(immune, "xgb", ncores) 
saveRDS(output, filename)
rm(immune)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="delta_xgb_third.rds"
delta <- dat$delta
output <- main(delta, "xgb", ncores) 
saveRDS(output, filename)
rm(delta)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="stellates_xgb_third.rds"
stellates <- dat$stellates
output <- main(stellates, "xgb", ncores) 
saveRDS(output, filename)
rm(stellates)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="endothelial_xgb_third.rds"
endothelial <- dat$endothelial
output <- main(endothelial, "xgb", ncores) 
saveRDS(output, filename)
rm(endothelial)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

##############################################################################################
