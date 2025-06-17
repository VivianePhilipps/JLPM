#' Fisher information function, used in Step 4 _ Selecting from the 4Smethod
#'
#' This function computes the contribution if each item during stages, based on 
#' the Fisher information. 
#' Specifically, this function xxx. 
#' The output is a list of two matrices : xxx.
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

IRT4Sselecting <- function(D,proj)
{
  if(missing(D)) stop("the model (argument D) is missing")
  if(!inherits(D,"jointLPM")) stop("use only with jointLPM model")
  if(missing(proj)) stop("the vector (argument proj) is missing")
  
  ns <- length(proj) + 1  # nb of stages
  ni <- D$N[12]  # nb of items
  nl <- max(D$nbmod)  # nb of modalities
    
  thres <- D$thres
  discrim <- D$discrim
  
  
  ### raw information
  raw_info <- matrix(NA, ncol = ns, nrow = ni)
  rownames(raw_info) <- D$Names$Ynames
  colnames(raw_info) <- as.roman(1:ns)
  
  for(i in 1:ni){ # item
    
    for(s in 1:ns){ # stage
      
      # boundaries of stage s
      if(s == 1){ # first stage
        proj.s.inf <- -Inf 
        proj.s.sup <- proj[1]
      }
      if(s == ns){ # last stage
        proj.s.inf <- proj[ns-1] 
        proj.s.sup <- +Inf
      }
      if(s != 1 && s != ns){ # other stages
        proj.s.inf <- proj[s-1] 
        proj.s.sup <- proj[s]
      }
      
      E1 <- 0 # first element = sum_i(integral)
      E2 <- 0 # second element = sum_i(substraction)
      for(m in 0:(nl-1)){ # modality
        E1 <- E1 + integrate(f=.integrate_info,
                             lower=proj.s.inf, upper=proj.s.sup,
                             i, thres, discrim, m, nl)$value
        E2 <- E2 + .Pim_prime(lambda=proj.s.sup,i,thres,discrim,m,nl) - .Pim_prime(lambda=proj.s.inf,i,thres,discrim,m,nl)
      }  
      
      raw_info[i,s] <- E1 - E2
      
    }
    
  }
  
  ### % of information
  percent_info <- raw_info
  
  for(i in 1:ni)
    percent_info[i,] <- percent_info[i,] / colSums(raw_info)
  percent_info <- percent_info * 100
  
  
  
  res <- list(raw_info=raw_info,percent_info=percent_info)
  
  return(res)
  
}
  