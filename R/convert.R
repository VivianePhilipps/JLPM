#' Conversion
#'
#' This function converts a jointLPM object to an object of type \code{multlcmm} or
#' \code{Jointlcmm} from lcmm package.
#' The conversion to \code{multlcmm} unable the use of the dedicated postfit and
#' predictions functions implemented in the lcmm package (plot, predictL, predictY,
#' predictlink, predictYcond, fitY, VarExpl).
#' The conversion to \code{Jointlcmm} permits the use of functions cuminc and
#' plot (with options "baselinerisk" or "survival").
#'
#' @param object an object of class \code{jointLPM}
#' @param to character. Either "multlcmm" or "Jointlcmm", indicating to which type the object should be converted.
#' @return an object of class \code{multlcmm} or \code{Jointlcmm}.
convert <- function(object, to)
{
    if(!(to %in% c("multlcmm", "Jointlcmm"))) stop("only conversions to multlcmm or Jointlcmm are implemented")

    z <- object$call
    z$maxiter <- 0
    z$verbose <- FALSE
    z$sharedtype <- NULL

    nrisqtot <- object$N[2]
    nvarxevt <- object$N[3]
    nasso <- object$N[4]
    nef <- object$N[5]
    ncontr <- object$N[6]
    nvc <- object$N[7]
    ncor <- object$N[9]
    ntrtot <- object$N[10]
    nalea <- object$N[11]
    ny <- object$N[12]
    
    if(to=="multlcmm")
    {
        z[[1]] <- as.name("multlcmm")
        z$survival <- NULL
        z$hazard <- NULL
        z$hazardrange <- NULL
        z$hazardnodes <- NULL
        z$TimeDepVar <- NULL
        z$logscale <- NULL
        z$startWeibull <- NULL

        index <- c(nrisqtot+nvarxevt+nasso + 1:(nef+ncontr+nvc+ncor),
                   nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea + 1:ny)
        if(nalea>0) index <- c(index, nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot + 1:nalea)
        index <- c(index, nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor + 1:ntrtot)

        B <- object$best[index]
    }
    else
    {
        z[[1]] <- as.name("Jointlcmm")
        z$methInteg <- NULL
        z$nMC <- NULL
        z$fixed <- formula(paste(object$Names$Ynames[1],"~",object$form$fixed[2]))
        z$randomY <- NULL
        if(length(object$call$startWeibull))
        {
           if(object$call$startWeibull != 0) warning("startWeibull is not taking into account")
        }
        if(object$linktype[1]==3)
        {
            warning("the outcome is treated as Gaussian")
            z$link <- "linear"
        }
        if(object$linktype[1]==0) z$link <- "linear"
        
        index <- c(1:(nrisqtot+nvarxevt),
                   nrisqtot+nvarxevt+nasso + 1:nef,
                   nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea + 1,
                   nrisqtot+nvarxevt+nasso+nef+ncontr + 1:(nvc+ncor))
        if(object$linktype[1]==3)
        {
            index <- c(index, nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor + 1:2)
        }
        else
        {
            index <- c(index, nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor + 1:ntrtot)
        }

        ## attention : dans jointLPM Var(u0)=1 alors que dans Jointlcmm Var(E)=1
        ## donc on divise la variance des EA par Var(E). Sinon pb avec chol(varcov).
        if(nvc>0) object$best[nrisqtot+nvarxevt+nasso+nef+ncontr + 1:nvc] <- object$best[nrisqtot+nvarxevt+nasso+nef+ncontr + 1:nvc]/object$best[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea + 1]^2
        object$best[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea + 1] <- 1/object$best[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea + 1]^2 
        
        B <- object$best[index]
    }

    z$B <- B
    m <- eval(z)
    m$conv <- object$conv
    
    varcov <- matrix(0, length(object$best), length(object$best))
    varcov[upper.tri(varcov, diag=TRUE)] <- object$V
    varcov <- t(varcov)
    varcov[upper.tri(varcov, diag=TRUE)] <- object$V
    v <- varcov[index,index]
    m$V <- v[upper.tri(v, diag=TRUE)]

    return(m)
}
