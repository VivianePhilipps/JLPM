#' @export
predictRE <- function(x, data, control = NULL, variance = FALSE)
{
    data <- as.data.frame(data)

    ## loglik arguments
    Fargs <- createFargs(x, data)

    ## no cor, no randomY
    if(Fargs$ncor0 > 0) stop("No implemented yet with cor")
    if(Fargs$nalea0 > 0) stop("No implemented yet with randomY")

    ## prepare result
    ns <- Fargs$ns0
    nea <- sum(Fargs$idea0)
    predRE <- matrix(NA, ns, 1 + nea)
    colnames(predRE) <- c(x$Names$ID, x$Names$Xnames[which(x$idea == 1)])
    varRE <- rep(NA, ns * nea * (nea + 1) / 2)

    ## id in (clean) data
    id <- as.numeric(Fargs$ID)
    Fargs <- Fargs[setdiff(names(Fargs), "ID")]

    ## remove unused arguments
    Fargs$ns0 <- NULL
    Fargs$methInteg0 <- NULL
    Fargs$nMC0 <- NULL
    Fargs$dimMC0 <- NULL
    Fargs$seqMC0 <- NULL
    Fargs$expectancy0 <- NULL
    if(Fargs$nfix0 > 0) # we don't need to distingish b from bfix
    {
        Fargs$npm0 <- Fargs$npm0 + Fargs$nfix0
        btot <- rep(NA, Fargs$npm0)
        btot[which(Fargs$fix0 == 0)] <- Fargs$b0
        btot[which(Fargs$fix0 == 1)] <- Fargs$bfix0
        Fargs$b0 <- btot
    }
    Fargs$fix0 <- NULL
    Fargs$nfix0 <- NULL
    Fargs$bfix0 <- NULL

    args_i <- Fargs

    ## optimization options
    defaultmla <- list(b = rep(0, nea), gr = NULL, hess = NULL, maxiter = 30, epsa = 0.0001, epsb = 0.0001, epsd = 0.0001, partialH = NULL, digits = 8, print.info = FALSE, blinding = TRUE, multipleTry = 25, nproc = 1, clustertype = NULL, file = "", .packages = NULL)
    controlmla <- control[which(names(control) %in% names(defaultmla))]
    defaultmla[names(controlmla)] <- controlmla
    
    args_mla <- c(list(m = nea, fn = REconddensity, minimize = FALSE),
                  defaultmla)

    ## different initial values for each subject
    matinit <- FALSE
    if(!is.null(controlmla$b))
    {
        binit <- controlmla$b
        if(is.matrix(binit))
        {
            if(nrow(binit) != ns) stop("control$b should be a vector or a matrix with as many rows as subjects")
            if(ncol(binit) != nea) stop("control$b should be a vector or a matrix with as many columns as random effects")

            matinit <- TRUE
        }
    }
    
    nmes <- Fargs$nmes0
#browser()
    ## loop over the subjects
    sumnmes <- 0
    sumvarRE <- 0
    for(i in 1:ns)
    {
        ## data of subject i
        args_i$Y0 <- Fargs$Y0[sumnmes + 1:sum(nmes[i, ])]
        args_i$X0 <- Fargs$X0[sumnmes + 1:sum(nmes[i, ]), ]
        args_i$Tentr0 <- Fargs$Tentr0[i]
        args_i$Tevt0 <- Fargs$Tevt0[i]
        args_i$Devt0 <- Fargs$Devt0[i]
        args_i$ind_survint0 <- Fargs$ind_survint0[i]
        args_i$nmes0 <- Fargs$nmes0[i, ]
        args_i$indiceY0 <- Fargs$indiceY0[sumnmes + 1:sum(nmes[i, ])]

        if(Fargs$idst0 %in% c(2, 4))
        {
            args_i$Xcl_Ti0 <- Fargs$Xcl_Ti0[i, ]
            args_i$Xcl_GK0 <- Fargs$Xcl_GK0[15 * (i - 1) + 1:15, ]
        }
        
        if(Fargs$idst0 %in% c(3, 4))
        {
            args_i$Xcs_Ti0 <- Fargs$Xcs_Ti0[i, ]
            args_i$Xcs_GK0 <- Fargs$Xcs_GK0[15 * (i - 1) + 1:15, ]
        }

        if(matinit)
            args_mla$b <- binit[i, ]

        ## find the mode of the distribution
        res_i <- do.call(marqLevAlg::mla, c(args_mla, args_i))

        ## save results
        if(res_i$istop %in% c(1, 3))
            predRE[i, ] <- c(id[i], res_i$b)
        else
            predRE[i, 1] <- id[i]

        ## variance
        if(variance == TRUE)
        {
            varRE[sumvarRE + 1:(nea * (nea + 1) / 2)] <- res_i$v
        }

        sumnmes <- sumnmes + sum(nmes[i, ])
        sumvarRE <- sumvarRE + nea * (nea + 1) / 2
    }

    if(variance == TRUE)
        return(list(predRE = predRE, varRE = varRE))
    else
        return(predRE)
}

#' @export
REconddensity <- function(ui0,Y0,X0,Tentr0,Tevt0,Devt0,ind_survint0,idea0,idg0,idcor0,idcontr0,
                   idsurv0,idtdv0,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0,ny0,
                   nv0,nobs0,nmes0,idiag0,ncor0,nalea0,npm0,b0,
                   epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,
                   idst0,nXcl0,Xcl_Ti0,Xcl_GK0,
                   Xcs_Ti0,Xcs_GK0,nonlin0,centerpoly0)
{
    res <- 0
    .Fortran(C_reconddensity,as.double(ui0),as.double(Y0),as.double(X0),as.double(Tentr0),as.double(Tevt0),as.integer(Devt0),as.integer(ind_survint0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(idcontr0),as.integer(idsurv0),as.integer(idtdv0),as.integer(typrisq0),as.integer(nz0),as.double(zi0),as.integer(nbevt0),as.integer(idtrunc0),as.integer(logspecif0),as.integer(ny0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(idiag0),as.integer(ncor0),as.integer(nalea0),as.integer(npm0),as.double(b0),as.double(epsY0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPLORD0),as.integer(idst0),as.integer(nXcl0),as.double(Xcl_Ti0),as.double(Xcl_GK0),as.double(Xcs_Ti0),as.double(Xcs_GK0),as.integer(nonlin0),as.double(centerpoly0),loglik_res=as.double(res))$loglik_res
}

#'@export
createFargs <- function(x, data)
{
    namesargs <- c("b0", "Y0", "X0", "Tentr0", "Tevt0", "Devt0", "ind_survint0",
                   "idea0", "idg0", "idcor0", "idcontr0", "idsurv0", "idtdv0",
                   "typrisq0", "nz0", "zi0", "nbevt0", "idtrunc0", "logspecif0",
                   "ny0", "ns0", "nv0", "nobs0", "nmes0", "idiag0", "ncor0",
                   "nalea0", "npm0", "nfix0", "bfix0", "epsY0", "idlink0",
                   "nbzitr0", "zitr0", "uniqueY0", "indiceY0", "nvalSPLORD0",
                   "fix0", "methInteg0", "nMC0", "dimMC0", "seqMC0", "idst0",
                   "nXcl0", "Xcl_Ti0", "Xcl_GK0", "Xcs_Ti0", "Xcs_GK0",
                   "nonlin0", "centerpoly0", "expectancy0")

    res <- vector("list", 1 + length(namesargs))
    names(res) <- c("ID", namesargs)

    best <- coef(x)
    posfix <- x$posfix
    fix <- rep(0, length(best))
    if(length(posfix))
    {
        res$b0 <- best[setdiff(1:length(best), posfix)]
        res$bfix0 <- best[posfix]
        res$nfix0 <- length(posfix)
        fix[posfix] <- 1
        
    }
    else
    {
        res$b0 <- best
        res$bfix0 <- 0
        res$nfix0 <- 0
    }
    res$fix0 <- fix
    res$npm0 <- length(res$b0)

    res$idea0 <- x$idea
    res$idg0 <- x$idg
    res$idcor0 <- x$idcor
    res$idcontr0 <- x$idcontr
    res$idsurv0 <- x$idsurv
    res$idtdv0 <- x$idtdv
    res$typrisq0 <- x$typrisq
    res$nz0 <- x$nz
    res$zi0 <- x$hazardnodes
    res$nbevt0 <- length(x$nevent)
    res$idtrunc0 <- ifelse(length(x$call$survival[[2]]) == 4, 1, 0)
    res$logspecif0 <- x$logspecif
    res$ny0 <- x$N[12]
    res$nv0 <- length(x$Names$Xnames)
    res$idiag0 <- x$idiag
    res$ncor0 <- x$N[9]
    res$nalea0 <- x$N[11]
    res$epsY0 <- x$epsY
    res$idlink0 <- x$linktype

    nbzitr <- rep(2, res$ny0)
    nbzitr[which(res$idlink0 == 2)] <- x$nbnodes
    res$nbzitr0 <- nbzitr
    
    res$zitr0 <- x$linknodes

    methInteg <- 3
    if(!is.null(x$call$methInteg))
    {
        if(as.character(x$call$methInteg) == "MCO") methInteg <- 1
        if(as.character(x$call$methInteg) == "MCA") methInteg <- 2
    }
    res$methInteg0 <- methInteg
    
    res$nMC0 <- x$nMC
    res$idst0 <- x$sharedtype
    res$nonlin0 <- x$nonlin
    res$centerpoly0 <- x$centerpoly
    res$expectancy0 <- 0

    if(!is.null(data))
    {
        resNA <- removeNA(x$terms, data)
        dataWithoutNA <- resNA$newdata

        nobsparY <- resNA$nmes
        uniqueY0 <- NULL
        indiceY0 <- NULL
        nvalSPLORD <- rep(0, res$ny0)
        sumnobsparY <- 0
        nb <- 0
        for(k in 1:res$ny0)
        {
            if((res$idlink0 != 2) & (res$idlink0 != 3))
            {
                indiceY0 <- c(indiceY0, rep(0, nobsparY[k]))
                next
            }

            yk <- dataWithoutNA[sumnobsparY + 1:nobsparY[k],]
            uniqueTemp <- sort(unique(yk))
            permut <- order(order(yk))

            if(length(as.vector(table(yk))) == length(uniqueTemp))
            {
                indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
                if(res$idlink0[k] == 2)
                {
                    indiceTemp <- nb + indice[permut]
                }
                else
                {
                    indiceTemp <- indice[permut]
                }
                
                nb <- nb + length(uniqueTemp)
                
                uniqueY0 <- c(uniqueY0, uniqueTemp)
                indiceY0 <- c(indiceY0, indiceTemp)
                nvalSPLORD[k] <- length(uniqueTemp)
            }
            else
            {
                uniqueY0 <- c(uniqueY0, yk)
                indiceY0 <- c(indiceY0, ifelse(res$idlink0[k] == 2, nb, 0) + c(1:length(yk)))
                nb <- nb + length(yk)
                nvalSPLORD[k] <- length(yk)
            }

            sumnobsparY <- sumnobsparY + nobsparY[k]
        }

        if(is.null(uniqueY0)) uniqueY0 <- 0        
        res$uniqueY0 <- uniqueY0
        res$indiceY0 <- indiceY0[order(dataWithoutNA[, x$Names$ID])]
        dataWithoutNA <- dataWithoutNA[order(dataWithoutNA[, x$Names$ID]),]
        res$Y0 <- as.vector(dataWithoutNA[, x$Names$Ynames])
        
        id <- unique(dataWithoutNA[, x$Names$ID])
        listnmes <- tapply(data[, x$Names$Ynames, drop = FALSE], data[, x$Names$ID], function(d) apply(d, 2, function(x) length(which(!is.na(x)))))
        if(length(x$Names$Ynames) > 1)
            res$nmes0 <- t(simplify2array(listnmes))
        else
            res$nmes0 <- matrix(listnmes, ncol = 1)
        
        formulas <- list(x$terms$fixed.right,
                         x$terms$random,
                         x$terms$surv.commun,
                         x$terms$surv.cause,
                         x$terms$form.cor)
        
        resX0 <- createX0(formulas, dataWithoutNA)
        res$X0 <- resX0$X0
        colnames(res$X0)[which(colnames(res$X0) == "(Intercept)")] <- "intercept"
        if(any(colnames(res$X0) != x$Names$Xnames)) stop("X0 !")

        res$nobs0 <- nrow(res$X0)
        maxmes <- max(res$nmes0)
        res$ns0 <- length(id)
        res$ID <- id

        if(res$nbevt0 > 0)
        {
            dataSurv <- dataWithoutNA[, c(x$Names$ID, x$Names$Tnames, x$Names$TimeDepVar.name)]
            dataSurv <- dataSurv[!duplicated(dataSurv[, x$Names$ID]),]

            if(res$idtrunc0 == 1)
            {
                res$Tentr0 <- dataSurv[, 2]
                res$Tevt0 <- dataSurv[, 3]
                res$Devt0 <- dataSurv[, 4]
            }
            else
            {
                res$Tentr0 <- rep(0, res$ns0)
                res$Tevt0 <- dataSurv[, 2]
                res$Devt0 <- dataSurv[, 3]
            }
            if(is.null(x$Names$TimeDepVar.name))
            {
                res$ind_survint0 <- rep(0, res$ns0)
            }
            else
            {
                Tint <- dataSurv[, x$Names$TimeDepVar.name]
                res$ind_survint0 <- (Tint < res$Tevt0) + 0
            }
            
            deriv <- ifelse(res$idst0 > 2, TRUE, FALSE)
            data_tmp <- dataWithoutNA[!duplicated(dataWithoutNA[, x$Names$ID]), setdiff(x$Names$Xvar, as.character(x$call$var.time)), drop = FALSE]
            X_GK <- matrixGK(data_tmp, x$terms$fixed.right, x$terms$random, as.character(x$call$var.time), res$idtrunc0, res$Tentr0, res$Tevt0, deriv, res$zi0[1,1])

            res$Xcl_Ti0 <- X_GK$Xpred_Ti 
            res$Xcl_GK0 <- X_GK$Xpred_cl 
            res$Xcs_Ti0 <- X_GK$Xpredcs_Ti 
            res$Xcs_GK0 <- X_GK$Xpred_cs 
            res$nXcl0 <- c(ncol(res$Xcl_Ti0), ncol(res$Xcl_GK0))

            if(is.null(res$Xcl_GK0)) res$Xcl_GK0 <- 0
            if(is.null(res$Xcs_Ti0)) res$Xcs_Ti0 <- 0
            if(is.null(res$Xcs_GK0)) res$Xcs_GK0 <- 0
        }
        else
        {
            res$Tentr0 <- rep(0, res$ns0)
            res$Tevt0 <- rep(0, res$ns0)
            res$Devt0 <- rep(0, res$ns0)
            res$ind_survint0 <- rep(0, res$ns0)
            res$nXcl0 <- c(1, 1)
            res$Xcl_Ti0 <- rep(0, res$ns0)
            res$Xcl_GK0 <- rep(0, 15 * res$ns0)
            res$Xcs_Ti0 <- rep(0, res$ns0)
            res$Xcs_GK0 <- rep(0, 15 * res$ns0)
        }
    }
    else
    {
        maxmes <- 20
    }
    
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3)
    {
        dimMC <- sum(res$idea0) + res$nalea0
        if(res$ncor0 > 0) dimMC <- dimMC + maxmes # maxmes dep des donnees
        if((dimMC > 0) & (res$nMC0 > 0))
        {
            sequnif <- spacefillr::generate_sobol_owen_set(res$nMC0, dimMC)
            seqMC <- apply(sequnif, 2, qnorm)
        }
    }
    res$dimMC0 <- dimMC
    res$seqMC0 <- seqMC

    return(res)
}
