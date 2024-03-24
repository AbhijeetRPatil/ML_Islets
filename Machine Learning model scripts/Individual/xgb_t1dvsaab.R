library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library("Matrix")
library("Seurat")
library("parallel")
## current working directory 
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/")

## Load all funcitons
source("/Libraries/libs/data_gen_t1dvsaab.R")

# set random seed
set.seed(333)

##############################################################################################################################
##------------------------------------------------- main resampling methods  ---------------------------------------------##
##############################################################################################################################

freq.matrix = function(train_matrix, test_matrix, method, ncores, test.x, train.y, test.y)
{
  ## XGB
  if(method=="xgb")
  {
    best_param = list()
    best_seednumber = 1441
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
    print(best_logloss_index)
    set.seed(best_seednumber)
    bst_model <- xgb.train(data=train_matrix, params=best_param, nrounds=nround, nthread=ncores)

    # Predict hold-out test set
    test_pred <- predict(bst_model, newdata = test_matrix)
    names(test_pred) <- rownames(test.x)
    df <- data.frame("cell_barcodes"= names(test_pred), test_pred)
    rownames(df) <- NULL
    test_pred <- ifelse (test_pred > 0.5,1,0)
    df$test_pred_class <- test_pred
    
    ## confusion matrix
    conf_mat <- confusionMatrix (as.factor(test_pred), as.factor(test.y))
    ## accuracy
    acc <- conf_mat$overall[[1]]
    ## sensitivity
    sens <- conf_mat$byClass[[1]]
    ## specificity
    spec <- conf_mat$byClass[[2]]
    ## metrics
    perf <- c(acc, sens, spec)    
    ## variable importance 
    total_genes <- colnames(train_matrix)
    importance_matrix = xgb.importance(feature_names = total_genes, model = bst_model)
    selected_genes <- importance_matrix$Feature
  }
  
  return(list(selected_genes, perf, df))
}

resampling = function(data, y, method, ncores)
{
  ## resampling
  ############## Split data into training and testing ###########################
  s.idx <- sample(1:nrow(data),as.integer(nrow(data)*.7),replace = F)
  
  train.x <- data[s.idx,]; train.y <- y[s.idx]
  test.x <-  data[-s.idx,]; test.y <-  y[-s.idx]
  
  ## Prepare train matrix
  train_matrix <- xgb.DMatrix(data = train.x, label = train.y)
  ## Prepare test matrix
  test_matrix <- xgb.DMatrix(data = test.x, label = test.y)
  
  ## xgb
  if(method=="xgb")
  {
    obj = freq.matrix(train_matrix, test_matrix, method, ncores, test.x, train.y, test.y)
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
  mat0 = list()
  mat1 = matrix(0,nrow=iter,ncol=3)
  mat2 = list()
  for(i in 1:iter)  ###########################################
  {
    # resampling
    obj0 = resampling(data,y,method,ncores)
    print(paste0("Iteration ", i, " is completed"))
    mat0[[i]] = obj0[[1]]
    print(paste0("The output of Iteration ", i, " is ", obj0[[2]]))
    mat1[i,] = obj0[[2]]
    print(paste0("The cell barcodes are written ", i, " is ", head(obj0[[3]])))
    mat2[[i]] = obj0[[3]]
  }
  metric.vec = apply(mat1,2,mean)
  sd.vec = apply(mat1,2,sd)
  
  return(list(mat0, mat1, metric.vec, sd.vec, mat2))
}

main = function(data, y, method, iter, ncores)
{
  results = resampling_results(data, y, method, iter, ncores)
  genes.vec     = results[[1]]
  mat1 = results[[2]]
  metric.vec = results[[3]]
  sd.vec = results[[4]]
  cell.df = results[[5]]
  
  return(list(genes.vec, mat1, metric.vec, sd.vec, cell.df))
}

start_time <- Sys.time()
print(paste0("Job sarted at ", start_time))
dat <- data.gen()
iter=100
ncores=30
filename="xgb_t1dvsaab.rds"
output <- main(dat[[1]], dat[[2]], "xgb", iter, ncores) 
saveRDS(output, filename)
end_time <- Sys.time()
print(paste0("Job completed at ", end_time))
total_time <- end_time - start_time
total_time
print(paste0("Total time taken is ", end_time))
