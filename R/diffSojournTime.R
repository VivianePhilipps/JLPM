#' Difference of sojourn times
#'
#' This function computes the difference between two sojourn times.
#'
#' @param x1 a sojournTime object, called with \code{draws = TRUE} and
#' \code{returndraws = TRUE}, computed for a specified seed.
#' @param x2 a sojournTime object, called with \code{draws = TRUE} and
#' \code{returndraws = TRUE}, computed with the same seed as \code{x1}.
#' 
#' @return returns the median, the 2.5\% and 97.5\% quantiles, the mean,
#' the standard deviation of the difference between the two sojourn times
#' and the number of removed draws (eventually due to computational issues) 

#' @export
diffSojournTime <- function(x1, x2)
{
    if(!inherits(x1, "sojournTime") | !inherits(x2, "sojournTime")) stop("Use only with sojournTime objects")
    if(!inherits(x1, "draws") | !inherits(x2, "draws")) stop("Please use returndraws = TRUE in the sojournTime call")
    if(!length(attr(x1, "seed")) | !length(attr(x1, "seed"))) stop("Please specify the seed in the sojournTime call")
    if(attr(x1, "seed") != attr(x2, "seed")) stop("The seeds should be the same in both sojournTime objects")
    if((length(attr(x1, "items")) == 1) & (length(attr(x2, "items")) == 1) & any(attr(x1, "items") != attr(x2, "items"))) stop("the two sojournTime objects should specify the same item in maxState")
    if(!all(attr(x1, "items") %in% attr(x2, "items")) & !all(attr(x2, "items") %in% attr(x1, "items"))) stop("the set of maxState items of one sojournTime object should be contained in the set of maxState items of the other sojournTime object")
    
    res <- x1 - x2
    res_50 <- quantile(res, probs=0.5, na.rm=TRUE)
    res_2.5 <- quantile(res, probs=0.025, na.rm=TRUE)
    res_97.5 <- quantile(res, probs=0.975, na.rm=TRUE)
    
    res_mean <- mean(res, na.rm=TRUE)
    res_sd <- sd(res, na.rm=TRUE)
    nn <- length(which(is.na(res)))

    result <- c(res_50, res_2.5, res_97.5, res_mean, res_sd, nn)
    class(result) <- "diffSojournTime"

    return(result)
}
