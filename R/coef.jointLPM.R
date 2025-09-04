#' @export
coef.jointLPM <- function(object, ...)
{
    ## vector of parameters
    b <- object$best

    avt <- sum(object$N[1:6])
    nvc <- object$N[7]
    nea <- sum(object$idea)
    idiag <- object$idiag
    
    ## replace varcov with cholesky
    if(nvc > 0)
    {
        if(idiag == 0)
        {
            b[avt + 1:nvc] <- object$cholesky[-1]
        }
        else
        {
            id <- 1:nea
            indice <- c(id + id * (id - 1) / 2)
            b[avt + 1:nvc] <- object$cholesky[indice[-1]]
        }

        names(b) <- sub("varcov", "cholesky", names(b))
    }

    return(b)
}
