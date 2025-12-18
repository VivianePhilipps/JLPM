## Function : createX0
## Args :
##  - form : a list of formula
##  - data : a data frame containing all variables mentionned in form
## Returns : a list containing
##  - X0 : the pooled (with cbind) model matrix. The number of lines is the same as in 'data'.
##  - idform : a list of the same length as 'form' specifying which term of X0 appear in each formula.
##
## Caution 1 : 'data' should not contain any missing values. There is no check for that within the function.
## Caution 2 : if interaction terms are present, the names are sorted. For example x:age becomes age:x in X0.
## Caution 3 : intercept is not systematically in the first column. It will only if form[[1]] contains an intercept.

#'@export
createX0 <- function(form, data)
{   
    nform <- length(form)

    ## combine all model.matrix
    X0 <- NULL
    namesxk <- vector("list", nform)
    for(k in 1:nform)
    {
        xk <- model.matrix(delete.response(terms(form[[k]])), data=data)
        if(ncol(xk) == 0) next;

        z <- strsplit(colnames(xk), split=":", fixed=TRUE)
        z <- lapply(z, sort)
        colnames(xk) <- sapply(z, function(x) paste(x, collapse=":"))
        namesxk[[k]] <- colnames(xk)
        
        X0 <- cbind(X0, xk)
    }
 
    ## remove duplicated columns
    if(!is.null(X0))
    {
        unames <- unique(colnames(X0))
        X0 <- X0[,unames, drop=FALSE]
        z.X0 <- strsplit(unames, split=":", fixed=TRUE)
        z.X0 <- lapply(z.X0, sort)
        namesX0 <- sapply(z.X0, function(x) paste(x, collapse=":"))
    }
    else
    {
        X0 <- matrix(0, nrow(data), 0)
    }
    
    ## indicator of presence in each formula
    idform <- vector("list", nform)
    for(k in 1:nform)
    {
        if(!is.null(namesxk[[k]]))
        {
            idform[[k]] <- (namesX0 %in% namesxk[[k]]) + 0
        }
        else
        {
            idform[[k]] <- rep(0, ncol(X0))
        }
    }

    return(list(X0=X0, idform=idform))
}
