####### Functions
accuracy <- function(truth, predicted)
{
  if(length(truth) > 0)
    sum(truth==predicted)/length(truth) else
      return(0)
}

sensitivity <- function(truth, predicted)
{
  if(sum(truth==1) > 0)
    sum(predicted[truth==1]==1)/sum(truth==1)   else
      return(0)
}

specificity <- function(truth, predicted)
{
  if(sum(truth==0) > 0)
    sum(predicted[truth==0]==0)/sum(truth==0)   else
      return(0)
}