library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library("Matrix")
library("Seurat")
library("parallel")
# library("pROC")
## current working directory 
# setwd("C:/Users/abhij/Desktop/ML_New/")
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/")
memory.limit(size=10000000)
## Load all functions
# source("C:/Users/abhij/Desktop/ML_New/celltypes_data_gen_t1dvsaab.R")
source("/Libraries/celltypes_data_gen_t1dvsaab_rna.R")
# set random seed
set.seed(666)
##############################################################################################################################
##------------------------------------------------- main resampling methods  ---------------------------------------------##
##############################################################################################################################

resampling = function(data, y, method,ncores)
{
  ex_df <- as.data.frame(y)
  ex_df <- cbind("cells"=rownames(ex_df), ex_df)
  ex_df$cells <- paste0(ex_df$`y`, "_", sub(".*-","",ex_df$cells))
  
  # Create empty vector to store predicted probabilities
  preds.resp <- list()
  preds.class <- list()
  accuracy <- list()
  test_pred <- list()
  out <- list()
  
  for (i in names(table(ex_df$cells))) 
  {
    print(i)
    train_y <- ex_df[ex_df$cells!=i,]
    test_y <- ex_df[ex_df$cells==i,]
    train_x <- data[rownames(train_y),]
    test_x <- data[rownames(test_y),]
    train_y$cells <- NULL
    test_y$cells <- NULL
    ## Prepare train matrix
    train_matrix <- xgb.DMatrix(data = train_x, label = train_y$y)
    ## Prepare test matrix
    test_matrix <- xgb.DMatrix(data = test_x, label = test_y$y)

    best_param = list()
    best_seednumber = 6666
    best_logloss = Inf
    best_logloss_index = 0
    
    for (iter in 1:10) 
    {
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
    print(best_logloss_index)
    set.seed(best_seednumber)
    bst_model <- xgb.train(data=train_matrix, params=best_param, nrounds=nround, nthread=ncores)
    
    preds.resp[[i]] <- predict(bst_model, as.matrix(test_x))
    preds.class[[i]] <- ifelse(preds.resp[[i]] > 0.5,1,0)
    accuracy[[i]] <- round(sum(preds.class[[i]] == test_y$y) / length(test_y$y)*100,2)
    
    # params <- list(
    #   objective = "binary:logistic",
    #   max_depth = 3,
    #   eta = 0.1,
    #   gamma = 0.1,
    #   subsample = 0.8,
    #   colsample_bytree = 0.8,
    #   eval_metric = "auc"
    # )
    # 
    # # Train XGBoost classifier on training data
    # xgb_model <- xgboost(
    #   data = train_matrix,
    #   # label = train_y,
    #   params = params,
    #   nrounds = 100,
    #   nthread = ncores,
    # )
    # 
    # preds.resp[[i]] <- predict(xgb_model, as.matrix(test_x)) 
    # preds.class[[i]] <- ifelse(preds.resp[[i]] > 0.5, 1, 0)
    # accuracy[[i]] <- round(sum(preds.class[[i]] == test_y$y) / length(test_y$y)*100,2)
  }
  return(list(preds.resp, preds.class, accuracy))
}

#################################################################
resampling_results = function(data, y, method, iter, ncores)
{
  mat0 = list()
  mat1 = matrix(0,nrow=iter,ncol=1)
  for(i in 1:iter)  ###########################################
  {
    # resampling
    obj0 = resampling(data, y, method, ncores)
    print(paste0("Iteration ", i, " is completed"))
  }
  return(list(obj0[[1]], obj0[[2]], obj0[[3]])) 
}

main = function(data,y,method,iter, ncores)
{
  results = resampling_results(data, y, method, iter, ncores)
  mu.vec     = results[[1]]
  mat1 = results[[2]]
  metric.vec = results[[3]]
  # sd.vec = results[[4]]
  return(list(mu.vec, mat1, metric.vec))
}
dat <- data.gen()
iter=1
ncores=30

## alpha
#start_time <- Sys.time()
#print(paste0("Job sarted at ", start_time))
#filename="alpha_xgb_t1dvsaab_03162024_cv_rna.rds"
#alpha <- dat$alpha
#output <- main(alpha[[1]], alpha[[2]], "xgb", iter, ncores)
#saveRDS(output, filename)
#rm(alpha)
#end_time <- Sys.time()
#print(paste0("Job completed at ", end_time))
#total_time <- end_time - start_time
#total_time
#print(paste0("Total time taken is ", end_time))

## beta
#start_time <- Sys.time()
#print(paste0("Job sarted at ", start_time))
#filename="beta_xgb_t1dvsaab_03162024_cv_rna.rds"
#beta <- dat$beta
#output <- main(beta[[1]], beta[[2]], "xgb", iter, ncores)
#saveRDS(output, filename)
#rm(beta)
#end_time <- Sys.time()
#print(paste0("Job completed at ", end_time))
#total_time <- end_time - start_time
#total_time
#print(paste0("Total time taken is ", end_time))

## acinar
#start_time <- Sys.time()
#print(paste0("Job sarted at ", start_time))
#filename="acinar_xgb_t1dvsaab_03162024_cv_rna.rds"
#acinar <- dat$acinar
#output <- main(acinar[[1]], acinar[[2]], "xgb", iter, ncores)
#saveRDS(output, filename)
#rm(acinar)
#end_time <- Sys.time()
#print(paste0("Job completed at ", end_time))
#total_time <- end_time - start_time
#total_time
#print(paste0("Total time taken is ", end_time))

## ductal
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="ductal_xgb_t1dvsaab_03162024_cv_rna.rds"
ductal <- dat$ductal
output <- main(ductal[[1]], ductal[[2]], "xgb", iter, ncores)
saveRDS(output, filename)
rm(ductal)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## immune
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="immune_xgb_t1dvsaab_03162024_cv_rna.rds"
immune <- dat$immune
output <- main(immune[[1]], immune[[2]], "xgb", iter, ncores) 
saveRDS(output, filename)
rm(immune)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

## delta
start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
filename="delta_xgb_t1dvsaab_03162024_cv_rna.rds"
delta <- dat$delta
output <- main(delta[[1]], delta[[2]], "xgb", iter, ncores) 
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
filename="endothelial_xgb_t1dvsaab_03162024_cv_rna.rds"
endothelial <- dat$endothelial
output <- main(endothelial[[1]], endothelial[[2]], "xgb", iter, ncores) 
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
filename="stellates_xgb_t1dvsaab_03162024_cv_rna.rds"
stellates <- dat$stellates
output <- main(stellates[[1]], stellates[[2]], "xgb", iter, ncores) 
saveRDS(output, filename)
rm(stellates)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))

