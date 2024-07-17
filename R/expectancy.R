## #' @export
## probaIRT <- function(model, t, Y, X, nmes=length(Y), indiceY, nMC=1000)
## {
##     ## calcule log(P(Y(t)<=y, T>t / X, theta))
    
##     ny <- model$N[12]
##     idlink <- model$linktype

##     Tentr <- 0
##     Devt <- 0
##     idea <- model$idea
##     idg <- model$idg
##     idcor <- model$idcor
##     idcontr <- model$idcontr
##     idsurv <- model$idsurv
##     idtdv <- model$idtdv
##     typrisq <- model$typrisq
##     nz <- model$nz
##     zi <- model$hazardnodes
##     nbevt <- length(model$nevent)
##     idtrunc <- 0
##     logspecif <- model$logspecif
##     nv <- length(model$Names$Xnames)
##     nobs <- length(Y)
##     nea <- sum(idea)
##     idiag <- model$idiag
##     ncor <- model$N[9]
##     nalea <- model$N[11]
##     epsY <- model$epsY
##     nbzitr <- rep(2,ny)
##     nbzitr[which(idlink==2)] <- model$nbnodes
##     zitr <- model$linknodes
##     fix <- rep(0,length(model$best))
##     posfix <- eval(model$call$posfix)
##     if(length(posfix)) fix[posfix] <- 1
##     methInteg <- 3
##     dimMC <- nea + nalea + sum(nmes)*as.numeric(ncor>0)
##     seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
##     npmtot <- length(model$best)
##     btot <- model$best

##     nvc <- model$N[7]
##     if(nvc>0)
##     {
##         ## remplacer varcov par cholesky dans btot

##         if(idiag==1)
##         {
##             btot[sum(model$N[1:6]) + 1:nvc] <- sqrt(btot[sum(model$N[1:6]) + 1:nvc])
##         }
##         else
##         {
##             btot[sum(model$N[1:6]) + 1:nvc] <- model$cholesky[-1]
##         }
##     }

    
##     ind_survint <- 0

##     nvalSPLORD <- rep(1,ny)
##     uniqueY <- NULL
##     for(k in 1:ny)
##     {
##         nvalSPLORD[k] <- length(model$mod[[k]])
        
##         if(idlink[k]==3) uniqueY <- c(uniqueY,model$mod[[k]])
##         if(idlink[k]==2) uniqueY <- c(uniqueY,model$zitr[1,k])        
##     }
    
        
##     proba <- 0

##     .Fortran(C_proba_irtsre,as.double(Y),as.double(X),as.double(Tentr),as.double(t),as.integer(Devt),as.integer(ind_survint),
##      as.integer(idea),as.integer(idg),as.integer(idcor),as.integer(idcontr),as.integer(idsurv),as.integer(idtdv),
##      as.integer(typrisq),as.integer(nz),as.double(zi),as.integer(nbevt),as.integer(idtrunc),as.integer(logspecif),
##      as.integer(ny),as.integer(nv),as.integer(nobs),as.integer(nea),as.integer(nmes),as.integer(idiag),as.integer(ncor),as.integer(nalea),
##      as.double(epsY),as.integer(idlink),as.integer(nbzitr),as.double(zitr),as.integer(uniqueY),as.integer(indiceY),
##      as.integer(nvalSPLORD),as.integer(fix),as.integer(methInteg),as.integer(nMC),as.integer(dimMC),as.double(seqMC),as.integer(npmtot),as.double(btot),res=as.double(proba))$res
## }


#' Computation of sojourn time from a joint model
#'
#' From a joint shared random effect model with ordinal longitudinal outcomes,
#' this function computes the time spend before reaching a certain level of the outcome(s).
#'
#' 1. Overall expected sojourn time
#'
#' The expected time before reaching level k from outcome Y is the integral over time t, from 0 to infinity, of P(Y(t) <= k, T > t), that is int_0^infty P(Y(t) <= k, T > t) dt.
#'
#' 2. Expected sojourn time from a time s
#'
#' Conditionally on being under level k at time s and being alive at time s, the sojourn time is int_s^infty P(Y(t) <= k, T > t | Y(s) <= k, T > s) dt = int_s^infty P(Y(t) <= k, T > t, Y(s) <= k) dt / P(Y(s) <= k, T > s)
#' 
#' 
#' @param x an object of class \code{jointLPM} representing of joint shared random effects model
#' @param event a named list specifying the state (ie the outcomes levels) in which the sojourn time is computed.
#' @param cond an optional named list specifying the initial state at start time (argument \code{start}).
#' @param newdata a data frame specifying the covariate profile from which the sojourn time is computed.
#' @param var.time a character string specifying the name of the time variable in the
#' longitudinal submodel. Note that this time covariate should not be included in newdata.
#' @param start a numeric value specifying the time from which the sojourn time is computed (the lower bound of the integral over time). Default to 0.
#' @param nMC an integer giving the number of Monte Carlo simulations used to compute the integral over the random effects.
#' @param upper a numeric specifying the upper bound of the integral over time. Default to 150.
#' @param subdivisions passed to the \code{integrate} function.
#' @param rel.tol passed to the \code{integrate} function.
#' @param draws logical indicating if confidence interval should be computed. Default to FALSE.
#' @param ndraws integer giving the number of draws used to compute the confidence interval. Default to 2000.
#' @param returndraws logical indicating if the \code{ndraws} results should be returned. Default to FALSE.
#' @param cl either a cluster created with \code{makeCluster} or an integer specifying the number of cores
#' that should be used for computation. Only used with draws = TRUE.
#'
#' @return if \code{draws = FALSE}, returns a single value. If \code{draws = TRUE} and \code{returndraws = FALSE},
#' returns the median, the 2.5\% and 97.5\% quantiles, the mean, the standard deviation and the number of removed draws (computational issues).
#' If \code{draws = TRUE} and \code{returndraws = TRUE}, returns the \code{ndraws} values.
#' 
#' @examples
#' library(lcmm)
#' paq <- paquid[which(paquid$age_init < paquid$agedem), ]
#' paq$age65 <- (paq$age - 65) / 10
#' paq$ageinit65 <- (paq$age_init - 65) / 10
#' paq$agedem65 <- (paq$agedem - 65) / 10
#' 
#' #### Estimation on a joint model with one ordinal marker
#' \dontrun{
#' M2 <- jointLPM(fixed = HIER ~ age65 * male,
#'                 random = ~ age65,
#'                 subject = "ID", 
#'                 link = "thresholds",
#'                 survival = Surv(ageinit65, agedem65, dem) ~ male,
#'                 sharedtype = 'RE',
#'                 var.time = "age65",
#'                 data = paq, 
#'                 methInteg = "QMC", 
#'                 nMC = 1000,
#'                 B = c(0.6, 2.399, -0.409, 2.076, 6.338, 0.994, -0.223, -0.005,
#'                      -0.299, 0.174, 0.523, 1.044, 1.064, 0.506))
#'
#' #### Computation of the expected time before reaching level 2 (HIER = 2)
#' #### (ie computes int_0^150 P(HIER(t) <= 2, T > t) dt)
#' expectancy(M2, list(HIER = 2), newdata = data.frame(male = 0), var.time = "age65")
#'
#' #### Computation of the expected time before reaching level 2 (HIER = 2),
#' #### when a level of at most 1 is reached at time 0.5
#' #### (ie computes int_0.5^150 P(HIER(t) <= 2, T > t | HIER(0.5) <= 1, T > 0.5) dt)
#' expectancy(M2, list(HIER = 2), cond = list(HIER = 1), start = 0.5,
#' newdata = data.frame(male = 0), var.time = "age65")
#' }
#' 
#' @export
#' 
expectancy <- function(x, event, cond=NULL, newdata, var.time, start=0, nMC=1000, upper=150, subdivisions=100L, rel.tol=.Machine$double.eps^0.25, draws=FALSE, ndraws=2000, returndraws=FALSE, cl=NULL)
{
    if(missing(x)) stop("the model (argument x) is missing")
    if(!inherits(x,"jointLPM")) stop("use only with jointLPM model")
    if(x$call$sharedtype == 'CL') stop("Not mplemented with current level (sharedtype = 'CL')")
    if(missing(event)) stop("argument event is missing")
#    if(!missing(cond) & (start==0)) stop("argument cond should only be used with start > 0")
    if(missing(var.time)) stop("argument var.time is missing")
    if(!(var.time %in% x$Names$Xvar)) stop("var.time does not appear in the model")
    if(missing(newdata) & length(setdiff(x$Names$Xvar,var.time))) stop("argument newdata is missing")
    if(!all(setdiff(x$Names$Xvar,var.time) %in% colnames(newdata))) stop(paste("newdata should include variables",paste(setdiff(x$Names$Xvar, var.time),collapse=" ")))

    if((x$conv != 1) & (draws != FALSE))
    {
        draws <- FALSE
        warning("The model did not converge properly. No confidence interval can be computed.")
    }
    
    if(any(x$idtdv==1))
    {
        if(is.na(newdata[,x$Names$TimeDepVar.name])) newdata[,x$Names$TimeDepVar.name] <- Inf
    }
    
    newdata1 <- na.omit(newdata[1,setdiff(x$Names$Xvar,c(var.time)),drop=FALSE])


    ny <- x$N[12]
    idlink <- x$linktype

    nmes <- rep(0,ny)
    nmescond <- rep(0,ny)
    Yevent <- rep(NA,ny)
    indiceYevent <- rep(NA,ny)
    Ycond <- rep(NA,ny)
    indiceYcond <- rep(NA,ny)
    nvalSPLORD <- rep(1,ny)
    uniqueY <- NULL
    indic_s1t2 <- rep(NA, 2*ny)
    for(k in 1:ny)
    {
        if(x$Names$Ynames[k] %in% names(cond))
        {
            if(idlink[k] != 3) stop("cond should only include ordinal outcomes")
            Ycond[k] <- cond[[x$Names$Ynames[k]]]
            nmes[k] <- nmes[k] + 1 # a faire : verifier qu'on a une seule valeur par Y
            nmescond[k] <- nmescond[k] + 1
            nvalSPLORD[k] <- length(x$mod[[k]])
            indiceYcond[k] <- which(x$mod[[k]] == Ycond[k])
            indic_s1t2[2*(k-1)+1] <- 1
        }
        
        if(x$Names$Ynames[k] %in% names(event))
        {
            if(idlink[k] != 3) stop("event should only include ordinal outcomes")
            Yevent[k] <- event[[x$Names$Ynames[k]]]
            nmes[k] <- nmes[k] + 1
            nvalSPLORD[k] <- length(x$mod[[k]])
            indiceYevent[k] <- which(x$mod[[k]] == Yevent[k])
            indic_s1t2[2*(k-1)+2] <- 2
        }

        if(idlink[k]==3) uniqueY <- c(uniqueY,x$mod[[k]])
        if(idlink[k]==2) uniqueY <- c(uniqueY,x$linknodes[1,k])
        
    }
    
    Y <- na.omit(as.vector(rbind(Ycond,Yevent)))
    indiceY <- na.omit(as.vector(rbind(indiceYcond,indiceYevent)))
    indic_s1t2 <- na.omit(indic_s1t2)
    
    fixed <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$fixed[2])
    random <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$random[2])
    contr <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$contr[2])
    surv <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.commun[2])
    survcause <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.cause[2])
    cor <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.cor[2])

    fixed <- formula(paste("~",fixed))
    random <- formula(paste("~",random))
    contr <- formula(paste("~",contr))
    surv <- formula(paste("~",surv))
    survcause <- formula(paste("~",survcause))
    cor <- formula(paste("~",cor))

    Tentr <- 0
    Devt <- 0
    idea <- x$idea
    idg <- x$idg
    idcor <- x$idcor
    idcontr <- x$idcontr
    idsurv <- x$idsurv
    idtdv <- x$idtdv
    typrisq <- x$typrisq
    nz <- x$nz
    zi <- x$hazardnodes
    nbevt <- length(x$nevent)
    idtrunc <- 0
    logspecif <- x$logspecif
    nv <- length(x$Names$Xnames)
    nobs <- length(Y)
    nea <- sum(idea)
    nvc <- x$N[7]
    idiag <- x$idiag
    ncor <- x$N[9]
    nalea <- x$N[11]
    epsY <- x$epsY
    nbzitr <- rep(2,ny)
    nbzitr[which(idlink==2)] <- x$nbnodes
    zitr <- x$linknodes
    fix <- rep(0,length(x$best))
    posfix <- eval(x$call$posfix)
    if(length(posfix)) fix[posfix] <- 1
    methInteg <- 3
    dimMC <- nea + nalea + sum(nmes)*as.numeric(ncor>0)
    seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
    npmtot <- length(x$best)
    btot <- x$best
    b <- btot[which(fix == 0)]
    bfix <- btot[which(fix == 1)]
    npm <- length(b)
    nfix <- length(bfix)

    if(nvc>0)
    {
        ## remplacer varcov par cholesky dans btot

        if(idiag==1)
        {
            btot[sum(x$N[1:6]) + 1:nvc] <- sqrt(btot[sum(x$N[1:6]) + 1:nvc])
        }
        else
        {
            btot[sum(x$N[1:6]) + 1:nvc] <- x$cholesky[-1]
        }
    }

    
    
    
    fctprob <- function(t, s, x, newdata, Y, fixed, random, contr, surv, survcause, cor,
                        indic_s1t2, Tentr,Devt,
                        idea,idg,idcor,idcontr,idsurv,idtdv,
                        typrisq,nz,zi,nbevt,idtrunc,logspecif,
                        ny,nv,nobs,nea,nmes,idiag,ncor,nalea,
                        epsY,idlink,nbzitr,zitr,uniqueY,indiceY,
                        nvalSPLORD,fix,methInteg,nMC,dimMC,seqMC,npm,b,nfix,bfix)
    {
        ## if(s>0)
        ## {
        ##     time <- unlist(lapply(as.vector(nmes),function(x) rep(c(s,t),length.out=x))) # suppose que si on a une seule mesure d'un Y, c'est au temps s
        ## }
        ## else
        ## {
        ##     time <- rep(t, sum(nmes))
        ## }
        time <- ifelse(indic_s1t2==1, s, t)
        ##cat("time=", time, "\n")
        newdata <- data.frame(newdata, t=time, row.names=NULL)
        
        Xfixed <- model.matrix(fixed, data=newdata)
        Xrandom <- model.matrix(random, data=newdata)
        Xcontr <- model.matrix(contr,data=newdata)
        Xsurv <- model.matrix(surv,data=newdata)
        Xsurvcause <- model.matrix(survcause,data=newdata)
        Xcor <- model.matrix(cor,data=newdata)   
        
        X0 <- cbind(Xfixed, Xrandom, Xsurv, Xsurvcause, Xcor)

        Xnames <- x$Names$Xnames
        Xnames[1] <- "(Intercept)"
        Xnames <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",Xnames)
        X0 <- X0[,Xnames,drop=FALSE]

        Tevt <- t
        ind_survint <- 0
        if(any(idtdv==1)) ind_survint <- as.numeric(newdata[1,x$Names$TimeDepVar.name] < t)

        expectancy <- 1
        proba <- 0
        
        res <- .Fortran(C_loglik,
                        as.double(Y),
                        as.double(X0),
                        as.double(Tentr),
                        as.double(Tevt),
                        as.integer(Devt),
                        as.integer(ind_survint),
                        as.integer(idea),
                        as.integer(idg),
                        as.integer(idcor),
                        as.integer(idcontr),
                        as.integer(idsurv),
                        as.integer(idtdv),
                        as.integer(typrisq),
                        as.integer(nz),
                        as.double(zi),
                        as.integer(nbevt),
                        as.integer(idtrunc),
                        as.integer(logspecif),
                        as.integer(ny),
                        as.integer(1),
                        as.integer(nv),
                        as.integer(nobs),
                        as.integer(nmes),
                        as.integer(idiag),
                        as.integer(ncor),
                        as.integer(nalea),
                        as.integer(npm),
                        as.double(b),
                        as.integer(nfix),
                        as.double(bfix),
                        as.double(epsY),
                        as.integer(idlink),
                        as.integer(nbzitr),
                        as.double(zitr),
                        as.double(uniqueY),
                        as.integer(indiceY),
                        as.integer(nvalSPLORD),
                        as.integer(fix),
                        as.integer(methInteg),
                        as.integer(nMC),
                        as.integer(dimMC),
                        as.double(seqMC),
                        as.integer(1),
                        as.integer(0),
                        as.double(0),
                        as.double(0),
                        as.integer(expectancy),
                        res=as.double(proba))$res

        if(!is.na(res)){ if(res == -1E-9) res <- NA }

        return(exp(res))
    }
    #browser()
    fctprobVect <- Vectorize(fctprob,"t")


    
    ## doone : computes the result for a set of parameters
    doone <- function(bdraw)
    {
        bdrawest <- bdraw[which(fix == 0)]
        bdrawfix <- bdraw[which(fix == 1)]
       ## b <- btot + Chol %*% bdraw
        res1 <- list(value=NA)
        try(res1 <- integrate(f=fctprobVect, lower=start, upper=upper, s=start, x=x,
                          newdata=newdata1, Y=Y,
                          fixed=fixed, random=random, contr=contr, surv=surv,
                          survcause=survcause, cor=cor,
                          indic_s1t2=indic_s1t2, Tentr=Tentr,Devt=Devt,
                          idea=idea,idg=idg,idcor=idcor,idcontr=idcontr,idsurv=idsurv,
                          idtdv=idtdv,
                          typrisq=typrisq,nz=nz,zi=zi,nbevt=nbevt,idtrunc=idtrunc,
                          logspecif=logspecif,
                          ny=ny,nv=nv,nobs=nobs,nea=nea,nmes=nmes,idiag=idiag,ncor=ncor,
                          nalea=nalea,
                          epsY=epsY,idlink=idlink,nbzitr=nbzitr,zitr=zitr,
                          uniqueY=uniqueY,indiceY=indiceY,
                          nvalSPLORD=nvalSPLORD,fix=fix,methInteg=methInteg,nMC=nMC,
                          dimMC=dimMC,seqMC=seqMC,npm=npm,b=bdrawest,nfix=nfix,bfix=bdrawfix,
                          rel.tol=rel.tol, subdivisions=subdivisions))

        if(class(res1)=="try-error") print(bdraw)
        
        result <- res1$value
  
        if((start > 0) & (length(na.omit(Ycond))))
        {
            nobscond <- length(na.omit(Ycond))
            Ycond <- na.omit(Ycond)
            indiceYcond <- na.omit(indiceYcond)
            
            res2 <- fctprob(t=start, s=0, x=x, newdata=newdata1, Y=Ycond,
                            fixed=fixed, random=random, contr=contr, surv=surv,
                            survcause=survcause, cor=cor,
                            indic_s1t2=rep(2, nobscond), Tentr=Tentr,Devt=Devt,
                            idea=idea,idg=idg,idcor=idcor,idcontr=idcontr,idsurv=idsurv,
                            idtdv=idtdv,
                            typrisq=typrisq,nz=nz,zi=zi,nbevt=nbevt,idtrunc=idtrunc,
                            logspecif=logspecif,
                            ny=ny,nv=nv,nobs=nobscond,nea=nea,nmes=nmescond,idiag=idiag,
                            ncor=ncor,nalea=nalea,
                            epsY=epsY,idlink=idlink,nbzitr=nbzitr,zitr=zitr,uniqueY=uniqueY,
                            indiceY=indiceYcond,
                            nvalSPLORD=nvalSPLORD,fix=fix,methInteg=methInteg,nMC=nMC,
                            dimMC=dimMC,seqMC=seqMC,
                            npm=npm,b=bdrawest,nfix=nfix,bfix=bdrawfix)
            
            result <- result/res2            
        }
 
        return(result)
    }


    
    if(!isTRUE(draws))
    {
        ## compute the result for parameters btot
        result <- doone(bdraw=btot)
    }
    else
    {
        ndraws <- as.integer(ndraws)
        
        posfix <- eval(x$call$posfix)
        
        if(ndraws>0)
        {
            Mat <- matrix(0,ncol=npmtot,nrow=npmtot)
            Mat[upper.tri(Mat,diag=TRUE)]<- x$V
            if(length(posfix))
            {
                Mat2 <- Mat[-posfix,-posfix]
                Chol2 <- chol(Mat2)
                Chol <- matrix(0,npmtot,npmtot)
                Chol[setdiff(1:npmtot,posfix),setdiff(1:npmtot,posfix)] <- Chol2
                Chol <- t(Chol)
            }
            else
            {
                Chol <- chol(Mat)
                Chol <- t(Chol)
            }
        }
        
        bdraw <- replicate(ndraws, btot + as.vector(Chol %*% rnorm(npmtot))) # ndraws colonnes

        ## compute the result for the ndraws sets of parameters
        if(!is.null(cl))
        {
            ncl <- NULL
            if(!inherits(cl,"cluster"))
            {
                if(!is.numeric(cl)) stop("argument cl should be either a cluster or a numeric value indicating the number of cores")
                
                ncl <- cl
                cl <- makeCluster(ncl)
            }
      
            ## set different seeds
            clusterSetRNGStream(cl)
            
            ## export other arguments
            clusterExport(cl, list("fctprobVect", "start", "upper", "x",
                          "newdata1", "Y",
                          "fixed", "random", "contr", "surv",
                          "survcause", "cor",
                          "indic_s1t2", "Tentr","Devt",
                          "idea", "idg", "idcor", "idcontr", "idsurv",
                          "idtdv",
                          "typrisq", "nz", "zi", "nbevt", "idtrunc",
                          "logspecif",
                          "ny", "nv", "nobs", "nea", "nmes", "idiag", "ncor",
                          "nalea",
                          "epsY", "idlink", "nbzitr", "zitr",
                          "uniqueY", "indiceY",
                          "nvalSPLORD", "fix", "methInteg", "nMC",
                          "dimMC", "seqMC", "npmtot", "btot",
                          "rel.tol", "subdivisions",
                          "fctprob", "Ycond", "indiceYcond"), envir = environment())
            
            ## get and export loaded packages
            pck <- .packages()
            dir0 <- find.package()
            dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
            clusterExport(cl,list("pck","dir"),envir=environment())
            clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))
            
            ## faire les ndraws repliques
            res <- parApply(cl,X=bdraw, MARGIN=2, FUN=doone)
            
            if(!is.null(ncl)) stopCluster(cl)            
        }
        else
        {
            res <- apply(bdraw, 2, doone)
        }
        
        if(returndraws) return(res)

        ## return the quantiles
        res_50 <- quantile(res, probs=0.5, na.rm=TRUE)
        res_2.5 <- quantile(res, probs=0.025, na.rm=TRUE)
        res_97.5 <- quantile(res, probs=0.975, na.rm=TRUE)

        res_mean <- mean(res, na.rm=TRUE)
        res_sd <- sd(res, na.rm=TRUE)
        nn <- length(which(is.na(res)))

        result <- c(res_50, res_2.5, res_97.5, res_mean, res_sd, nn)
    }

    return(result)
}



