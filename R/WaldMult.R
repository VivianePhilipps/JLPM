################################################################################
##                                                                            ##
##                               WaldMult                                     ##
##                                                                            ##
################################################################################

## Performs multivariate Wald tests
##
## Given a vector of estimated parameters (coef = (theta_1, ... theta_n)') and
## their associated variance matrix (vcov = Var(coef)), the WaldMult function
## performs either the test of the null hypothesis :
##
## 1. H0 : theta_m = (theta_j1, ..., theta_jm)' = 0
##    ie, tests if a subset of m parameters are simultaneously equal to zero.
##    The Wald statistic is w = theta_m' * Var(theta_m)^(-1) * theta_m
##    and the associated p-value is obtained by p = P(chi2_{dof=m} > w)
##
## 2. H0 : theta = c1*theta_j1 + ... + cm*theta_jm = 0
##    ie, tests if a combination of m parameters is equal to zero.
##    The Wald statistic is w = theta^2 / Var(theta)
##                            = (c1*theta_j1 + ... + cm*theta_jm)^2 / Var(c1*theta_j1 + ... + cm*theta_jm)
##    and the associated p-value is obtained by p = P(chi2_{dof=1} > w)
##
## 3. H0 : Theta = A*coef = 0 with A a matrix with m lines and n columns.
##    ie, tests if multiple combinations of parameters are simultaneously equal to zero.
##    The Wald statistic is w = Theta * Var(Theta)^(-1) * Theta'
##                            = (A*coef)' * (A*vcov*A')^(-1) * (A*coef)
##    and the associated p-value is obtained by p = P(chi2_{dof=m} > w)
##
##
## Parameters :
## coef : the vector of estimated parameters
## vcov : the matrix of variance of coef
## pos : vector containing the indices of the parameters involved in the test. Required if A is not specified.
## contrasts : optional vector containing the the c1, ..., cm value in case 2. Should be of same lengtht as pos.
## A : matrix A of case 3. A should have as many columns as parameters in coef.
## name : characters string giving the name of the test printed in the output (the row names of the output).
##        By default, the name's test is the null hypothesis.
## value : the value against which to test. By defualt, value=0.
##
## Returns :
## If contrasts is NULL, the function returns a matrix with 1 row and
## 2 columns containing the value of the Wald test's statistic and
## the associated p-value.
##
## If contrasts is not NULL, the function returns a matrix with 1 row
## and 4 columns containing the value of the coefficient (dot product
## of pos and contrasts), his standard deviation, the value of the
## Wald test's statistic and the associated p-value.
##
## If A is specified, the function returns a matrix with m rows and
## 2 columns containing the value of the Wald test's statistic and
## the associated p-value.
##
## Depends : none
##
## Author : Viviane Philipps
##
## Example :
##
## library(nlme) # for Orthodont data
## library(lcmm) # for hlme function
##
## Orthodont$ID <- as.numeric(Orthodont$Subject)
##
## m <- hlme(distance~age+Sex, data=Orthodont, random=~1, subject="ID")
## summarry(m)
##
## mm <- hlme(distance~-1+age+Sex, data=Orthodont, random=~1, subject="ID")
## summary(mm)
##
## Retrieve from model mm the difference between gender as in model m :
## WaldMult(coef(mm), vcov(mm), pos=c(2,3), contrasts=c(-1,1))
##  or
## WaldMult(coef(mm), vcov(mm), A=matrix(c(0,-1,1,0,0),1,5))
##
#' @export
WaldMult <- function(Mod, coef, vcov, pos = NULL, contrasts = NULL, A = NULL, name = NULL, value = 0)
{
    if(missing(Mod))
    {
        if(missing(coef)) stop("coef should be specified")
        if(missing(vcov)) stop("vcov should be specified")
        if(!is.matrix(vcov)) stop("vcov should be a matrix")
    }
    else
    {
##        if(!inherits(Mod, "jointLPM")) stop("Mod should be a jointLPM model")
        coef <- stats::coef(Mod)
        vcov <- stats::vcov(Mod)
    }
    
    npm <- length(coef)
    if((nrow(vcov) != npm) | (ncol(vcov) != npm)) stop("the dimension of vcov does not match with the length of coef")

    if(is.null(names(coef))) names(coef) <- paste("coef", 1:npm, sep = "")
    
    if(is.null(A)) ## cas 1. ou 2.
    {
        ## On teste si pos est un vecteur
        if (is.null(pos))
        {
            stop("Either A or pos must be specified")
        }
        else
        {
            if (!is.vector(pos)) stop("pos must be a numeric vector")

            ## subset of parameters and their variance
            Vect <- coef[pos]
            Mat <- vcov[pos, pos]
            
            
            if (is.null(contrasts)) ## cas 1.
            { 
                
                if (!missing(value))
                {
                    if (!is.vector(value)) stop("value must be a numeric vector")
                    if (length(value) != length(pos)) stop("value must have the same length as the vector pos")
                    
                    Vect <- coef[pos] - value
                }
                
                
                Wald <- t(Vect) %*% solve(Mat) %*% Vect
                
                ddl <- length(pos)
                p_value <- 1 - pchisq(Wald, df = ddl)
                
                Results <- matrix(NA, nrow = 1, ncol = 2)
                colnames(Results) <- c("Wald Test", "p_value")
                if (is.null(name)) 
                {
                    if (!is.null(value))
                    {
                        rownames(Results) <- paste(names(coef[pos]), " = ", value, collapse = " and ", sep = "")
                    }
                    else
                    {
                        rownames(Results) <- paste(paste(names(coef[pos]), collapse = " = "), "= 0")
                    }
                }
                else
                {
                    rownames(Results) <- name
                }
                
                Results[, 1] <- round(Wald, 5)
                Results[, 2] <- round(p_value, 5)
            }
            else ## cas 2.
            {
                if (length(contrasts) != length(pos))  stop("contrasts must have the same length as the vector pos")
                if (sum(abs(contrasts)) == 0) stop("The absolute value of the sum of contratsts components must be different from 0")
                
                ## linear combination of the parameters
                Scalaire <- sum(Vect * contrasts)
                
                if (value != 0)
                {
                    if (!is.vector(value)) stop("value must be a numeric vector")
                    if (length(value) != 1) stop("value must be a vector with a unique argument")
                    
                    Scalaire <- sum(Vect * contrasts) - value
                }

                ## variance of the combination
                Var <- t(contrasts) %*% Mat %*% contrasts
                

                Wald <- Scalaire^2 / Var
                p_value <- 1 - pchisq(Wald, df = 1) #2 * (1 - pnorm(abs(Wald)))
                
                Results <- matrix(NA, nrow = 1, ncol = 4)
                colnames(Results) <- c("coef", "Se", "Wald Test", "p_value")
                if (is.null(name)) 
                {
                    if(is.null(value)) value <- 0
                    rownames(Results) <- paste(paste(names(coef[pos]), "*", contrasts, collapse = " + "), "= ", value)
                } 
                else
                {
                    rownames(Results) <- name
                }
                
                
                Results[, 1] <- round(sum(Vect * contrasts), 5)
                Results[, 2] <- round(sqrt(Var), 5)
                Results[, 3] <- round(Wald, 5)
                Results[, 4] <- round(p_value, 5)
            }
            
        }
    }
    else ## cas 3.
    {
        if(ncol(A) != npm) stop("A should have as many columns as parameters in coef")
        ddl <- nrow(A)

        if(!missing(value))
        {
            if(length(value) == 1)
            {
                value <- rep(value, ddl)
            }
            else
            {
                if(length(value) != ddl) stop(paste("Vector value should be of length", ddl))
            }
        }
        else
        {
            value <- rep(0, ddl)
        }

        Theta <- A %*% coef - value
        Wald <- t(Theta) %*% solve(A %*% vcov %*% t(A)) %*% Theta
        p_value <- 1 - pchisq(Wald, df = ddl)
        
        Results <- matrix(NA, nrow = 1, ncol = 2)
        colnames(Results) <- c("Wald Test", "p_value")   
        if (is.null(name)) 
        {
            rownames(Results) <- paste("A*coef = (", paste(value, collapse = ", "), ")'", sep = "")
        }
        else
        {
            rownames(Results) <- name
        }            
        Results[, 1] <- round(Wald, 5)
        Results[, 2] <- round(p_value, 5)
    }
    
    return(Results)
}


