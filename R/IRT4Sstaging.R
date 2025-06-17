#' Projection function, used in Step 3 _ Staging from the 4Smethod
#'
#' This function determines the stage thresholds in the latent scale of the
#' dimension process. Specifically, this function predicts the dimension scores
#' corresponding to stage changes and deducts the corresponding thresholds in the
#' dimension scale based on a correspondence between the predicted scores and
#' the expected sum of the items. The output is a vector of the estimated stage transition
#' thresholds. 
#' 
#' @param D an object of class \code{jointLPM} representing the estimated 
#' dimension model with the items
#' @param Dss an object of class \code{jointLPM} representing the estimated proxy 
#' model between the dimension scores and the stages
#' @param nsim number of points used in the numerical integration (Monte-Carlo). 
#' nsim should be relatively important (nsim=1000 by default).
#' @param bounds dimension scale boundaries between which are computed a grid of
#' the dum of the items
#'
#' @return a vector of the estimated stage transition thresholds
#' 
#' @examples
#' 
#' @author Tiphaine Saulnier, Viviane Philipps and Cecile Proust-Lima
#' 
#' @references
#' preprint :  arXiv:2407.08278
#' 
#' @export

IRT4Sstaging <- function(D,Dss,nsim=1000,bounds=c(-10,10))
{
  if(missing(D)) stop("the model (argument D) is missing")
  if(!inherits(D,"jointLPM")) stop("use only with jointLPM model")
  if(missing(Dss)) stop("the model (argument Dss) is missing")
  if(!inherits(Dss,"jointLPM")) stop("use only with jointLPM model")
  
  
  # convert JLPM -> multlcmm for predictions
  dss <- convert(object = Dss, to = "multlcmm")
  
  
  ### item score predictions  
  Ypred <- predictYcond(x = dss,
                        lprocess = Dss$thres[,1],
                        nsim=nsim)$pred$Ypred[1:length(Dss$thres[,1])]
  
  ### conditional expectation of the sum of items
  tmp <- data.frame(lambda = seq(from = bounds[1], to = bounds[2], by = 0.01))
  tmp$sum <- sapply(X = tmp$lambda, 
                    FUN = .expect, 
                    ni=D$N[12], nl=max(D$nbmod), discrim = D$discrim, thres = D$thres)

  ### deduct lambda, such that score = sum of items
  res <- rep(NA,length(Dss$thres[,1]))
  for(j in 1:length(Dss$thres[,1]))
    res[j] <- max(tmp$lambda[which(abs(tmp$sum - Ypred[j]) == min(abs(tmp$sum - Ypred[j])))]) 
  
  
  return(res)
}
