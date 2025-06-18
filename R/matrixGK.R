matrixGK <- function(data, fixed, random = NULL, var.time, idtrunc, T0 = NULL, T, deriv = FALSE)
{
    ## ajouter variable temps aux donnees
    data$tmp_name <- NA
    colnames(data)[which(colnames(data) == "tmp_name")] <- var.time
    data[, var.time] <- T

    ## matrices pour risque au temps Ti
    mat_ef <- model.matrix(fixed, data = data)
    mat_ef <- as.data.frame(mat_ef[, -1])  # sans intercept car non-estime
    mat_ea <- NULL
    if(!is.null(random)) mat_ea <- model.matrix(random, data = data)
    Xpred_Ti <- as.matrix(cbind(data[, var.time], mat_ef, mat_ea))

    Xdt_Ti <- NULL
    if(deriv)
    {
        h <- sapply(data[, var.time], function(x) { max(1E-7, 1E-4 * abs(x)) })

        ## t + h
        dplus <- data
        dplus[, var.time] <- dplus[, var.time] + h
        mat_ef <- model.matrix(fixed, data = dplus)
        mat_ef <- as.data.frame(mat_ef[, -1])
        mat_ea <- NULL
        if(!is.null(random)) mat_ea <- model.matrix(random, data = dplus)
        Xplus <- cbind(mat_ef, mat_ea)

        ## t - h
        dmoins <- data
        dmoins[, var.time] <- dmoins[, var.time] - h
        mat_ef <- model.matrix(fixed, data = dmoins)
        mat_ef <- as.data.frame(mat_ef[, -1])
        mat_ea <- NULL
        if(!is.null(random)) mat_ea <- model.matrix(random, data = dmoins)
        Xmoins <- cbind(mat_ef, mat_ea)
        
        ## derivee par rapport a var.time
        Xdt_Ti <- as.matrix((Xplus - Xmoins) / (2 * h))
        Xdt_Ti <- derivMat(fixed, random, data, var.time)
    }

    ## points de quadrature Gauss-Kronrod
    ptGK_1 <- 0.991455371120812639206854697526329
    ptGK_2 <- 0.949107912342758524526189684047851
    ptGK_3 <- 0.864864423359769072789712788640926
    ptGK_4 <- 0.741531185599394439863864773280788
    ptGK_5 <- 0.586087235467691130294144838258730
    ptGK_6 <- 0.405845151377397166906606412076961
    ptGK_7 <- 0.207784955007898467600689403773245
    ptGK_8 <- 0
    
    ptsGK <- c(ptGK_1, -ptGK_1, ptGK_2, -ptGK_2, ptGK_3, -ptGK_3, ptGK_4, -ptGK_4, ptGK_5, -ptGK_5, ptGK_6, -ptGK_6, ptGK_7, -ptGK_7, ptGK_8) # integration [-1, 1]


    ## 15 lignes par sujet
    data$T0 <- T0
    data$T <- T
    data <- data[rep(1:nrow(data), each = 15), ] 

    ## matrices pour survie 0->T
    data[, var.time] <- data$T / 2 + data$T / 2 * ptsGK # integration [0, T]
    mat_ef <- model.matrix(fixed, data = data)
    mat_ef <- as.data.frame(mat_ef[,-1])
    mat_ea <- NULL
    if(!is.null(random)) mat_ea <- model.matrix(random, data = data)
    Xpred <- as.matrix(cbind(data[, var.time], mat_ef, mat_ea))

    Xpred0 <- NULL
    if(idtrunc)
    {
        ## matrices pour survie 0->T0
        data[, var.time] <- data$T0 / 2 + data$T0 / 2 * ptsGK # integration [0, T0]
        mat_ef <- model.matrix(fixed, data = data)
        mat_ef <- as.data.frame(mat_ef[, -1])
        mat_ea <- NULL
        if(!is.null(random)) mat_ea <- model.matrix(random, data = data)
        Xpred0 <- as.matrix(cbind(data[, var.time], mat_ef, mat_ea))
    }


    Xdt <- NULL
    Xdt0 <- NULL
    if(deriv)
    {
        data[, var.time] <- rep(T, each = 15)
        data[, var.time] <- data$T / 2 + data$T / 2 * ptsGK
        
        h <- sapply(data[, var.time], function(x) { max(1E-7, 1E-4 * abs(x)) })

        ## t + h
        dplus <- data
        dplus[, var.time] <- dplus[, var.time] + h
        mat_ef <- model.matrix(fixed, data = dplus)
        mat_ef <- as.data.frame(mat_ef[, -1])
        mat_ea <- NULL
        if(!is.null(random)) mat_ea <- model.matrix(random, data = dplus)
        Xplus <- cbind(mat_ef, mat_ea)

        ## t - h
        dmoins <- data
        dmoins[, var.time] <- dmoins[, var.time] - h
        mat_ef <- model.matrix(fixed, data = dmoins)
        mat_ef <- as.data.frame(mat_ef[, -1])
        mat_ea <- NULL
        if(!is.null(random)) mat_ea <- model.matrix(random, data = dmoins)
        Xmoins <- cbind(mat_ef, mat_ea)

        ## derivee par rapport a var.time
        Xdt <- as.matrix((Xplus - Xmoins) / (2 * h))

        Xdt <- derivMat(fixed, random, data, var.time)

        Xdt0 <- NULL
        if(idtrunc)
        {
            data[, var.time] <- rep(T0, each = 15)
            data[, var.time] <- data$T0 / 2 + data$T0 / 2 * ptsGK
            
            h <- sapply(data[, var.time], function(x) { max(1E-7, 1E-4 * abs(x)) })
            
            ## t + h
            dplus <- data
            dplus[, var.time] <- dplus[, var.time] + h
            mat_ef <- model.matrix(fixed, data = dplus)
            mat_ef <- as.data.frame(mat_ef[, -1])
            mat_ea <- NULL
            if(!is.null(random)) mat_ea <- model.matrix(random, data = dplus)
            Xplus <- cbind(mat_ef, mat_ea)
            
            ## t - h
            dmoins <- data
            dmoins[, var.time] <- dmoins[, var.time] - h
            mat_ef <- model.matrix(fixed, data = dmoins)
            mat_ef <- as.data.frame(mat_ef[, -1])
            mat_ea <- NULL
            if(!is.null(random)) mat_ea <- model.matrix(random, data = dmoins)
            Xmoins <- cbind(mat_ef, mat_ea)
            
            ## derivee par rapport a var.time
            Xdt0 <- as.matrix((Xplus - Xmoins) / (2 * h))
            Xdt0 <- derivMat(fixed, random, data, var.time)
        }
    }

    res <- list(Xpred_Ti = Xpred_Ti, Xpred_cl = cbind(Xpred, Xpred0), Xpredcs_Ti = Xdt_Ti, Xpred_cs = cbind(Xdt, Xdt0))
    
    return(res)
}
