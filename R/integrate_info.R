# internal function for the 4S method

#### STEP 4 - SELECTING
## function to integrate : (P'im)^2 / Pim between boundaries of stage s
.integrate_info <- function(lambda,
                              i,thres,discrim,
                              m,nl){
  # arguments : - lamba = dimension process value
  #             - i = item
  #             - thresh = item threshold parameters
  #             - discrim = item discrimination parameters
  #             - m = modality
  #             - nl = nb of levels
  res <- .Pim_prime(lambda,i,thres,discrim,m,nl)^2
  denom <- .Pim(lambda,i,thres,discrim,m,nl)
  res <- ifelse(denom == 0,
                res <- 0,
                res <- res / denom)
  return(res)
}
