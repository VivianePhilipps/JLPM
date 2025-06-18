derivMat <- function(fixed, random, data, var.time)
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
    Xdt <- as.matrix((Xplus - Xmoins) / (2 * h))

    return(Xdt)
}
