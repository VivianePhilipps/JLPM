# internal function for the 4S method

#### STEP 3 - STAGING
## conditional expectation function : E[sum(Y)] = sum_k{M - sum(P(Yk<=m))_m=0,1,...,m-1}
.expect <- function(lambda, ni, nl, discrim, thres){
  proba <- matrix(NA,nrow=nl-1,ncol=ni) # probability matrix P(Y<=m) with m=0,1,2,..,(nl-2)
  for(k in 1:ni)
    for(l in 1:(nl-1))
      proba[l,k] <- pnorm(-discrim[k]*(lambda-thres[l,k])) # gaussian distrib fct
  probaSum <- colSums(proba) # vector of sums for each item
  probaSum <- sum(probaSum)
  res <- ni*(nl-1) - probaSum
  return(res)
}
