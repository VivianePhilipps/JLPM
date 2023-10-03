#' Standard methods for estimated models
#' 
#' coef and vcov  methods for estimated jointLPM models.
#' 
#' 
#' @aliases coef.jointLPM vcov.jointLPM
#' @param object an object of class \code{jointLPM}
#' @param ...  other arguments. There are ignored in these functions.
#' @return For \code{coef}, the vector of the estimated parameters.
#' 
#' For \code{vcov}, the variance-covariance matrix of the estimated parameters.
#' 
#' @author Viviane Philipps
#' 
#' @name StandardMethods





#' @export
coef.jointLPM <- function(object, ...)
{
    
 if(object$conv %in% c(1,2,3))
 {
     res <- object$best

     if(object$N[7]>0)
     {
         res[sum(object$N[2:6])+1:object$N[7]] <- object$cholesky[-1]
         
         names(res) <- sub("varcov","cholesky",names(res))
     }
 }
 else
 {
  res <- NA
  cat("Output can not be produced since the program stopped abnormally. \n")
 }
 
 return(res)
}


#' @export
vcov.jointLPM <- function(object, ...)
{
    
 if(object$conv %in% c(1,2,3))
 {
  res <- matrix(0,length(object$best),length(object$best))
  res[upper.tri(res,diag=TRUE)] <- object$V
  res <- t(res)
  res[upper.tri(res,diag=TRUE)] <- object$V

  noms <- sub("varcov","cholesky",names(object$best))
  colnames(res) <- noms
  rownames(res) <- noms
 }
 else
 {
  res <- NA
  cat("Output can not be produced since the program stopped abnormally. \n")
 }

 return(res)
}
