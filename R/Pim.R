# internal function for the 4S method

#### STEP 4 - SELECTING
## probability function P : Pim(lambda) = P(Yi = m |lambda)
.Pim <- function(lambda,
                 i,thres,discrim,
                 m,nl){
  # arguments : - lamba = dimension process value
  #             - i = item
  #             - thresh = item threshold parameters
  #             - discrim = item discrimination parameters
  #             - m = modality
  #             - nl = nb of levels
  res <- NA
  discrim_i <- discrim[i]
  # m=0, lower modality
  if(m == 0){
    thres.im <- thres[1,i]
    res <- pnorm(-discrim_i * (lambda - thres.im))
  }
  # m=nl-1, upper modality
  if(m == (nl-1)){
    thres.im <- thres[nl-1,i]
    res <- 1 - pnorm(-discrim_i * (lambda - thres.im))
  }
  # other modalities
  if(m != 0 && m!=(nl-1)){
    thres.im.inf <- thres[m,i]
    thres.im.sup <- thres[m+1,i]
    res <- pnorm(-discrim_i * (lambda - thres.im.sup)) - pnorm(-discrim_i * (lambda - thres.im.inf))
  }
  return(res)  
}