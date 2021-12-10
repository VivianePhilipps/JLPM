#' @export
Mspline <- function(t,z)
{
    nz <- length(z)
    z <- sort(z)
    if(any(t<z[1]) | any(t>z[nz])) stop("t en dehors de z")
    
    ## repeter les noeuds externes
    zz <- rep(0,nz+6)
    zz[1:3] <- rep(z[1],3)
    zz[3+1:nz] <- z
    zz[3+nz+1:3] <- z[nz]
    
    ## declaration du resultat
    Tmm <- matrix(0,nrow=length(t),ncol=4)
    Tim <- matrix(0,nrow=length(t),ncol=4)
        
    ## boucle sur les temps
    for(j in 1:length(t))
    {
        ## encadrer t par les noeuds
        tz <- c(z[-c(1,nz)],t[j])
        itz <- c(rep(0,nz-2),1)
        l <- 3+which(itz[order(tz)]==1)
        
        ## calcul des Tmm et Tim non nuls
        
        ht <- t[j]-zz[l]
        htm <- t[j]-zz[l-1]
        h2t <- t[j]-zz[l+2]
        ht2 <- zz[l+1]-t[j]
        ht3 <- zz[l+3]-t[j]
        hht <- t[j]-zz[l-2]
        h <- zz[l+1]-zz[l]
        hh <- zz[l+1]-zz[l-1]
        h2 <- zz[l+2]-zz[l]
        h3 <- zz[l+3]-zz[l]
        h4 <- zz[l+4]-zz[l]
        h3m <- zz[l+3]-zz[l-1]
        h2n <- zz[l+2]-zz[l-1]
        hn <- zz[l+1]-zz[l-2]
        hh3 <- zz[l+1]-zz[l-3]
        hh2 <- zz[l+2]-zz[l-2]
        
        if(t[j]<z[nz])
        {           
            Tmm[j,4] <- (4*ht2*ht2*ht2)/(h*hh*hn*hh3)
            Tmm[j,3] <- (4*hht*ht2*ht2)/(hh2*hh*h*hn) + (-4*h2t*htm*ht2)/(hh2*h2n*hh*h) + (4*h2t*h2t*ht)/(hh2*h2*h*h2n)
            Tmm[j,2] <- 4*(htm*htm*ht2)/(h3m*h2n*hh*h) + (-4*htm*ht*h2t)/(h3m*h2*h*h2n) + (4*ht3*ht*ht)/(h3m*h3*h2*h)
            Tmm[j,1] <- 4*(ht*ht*ht)/(h4*h3*h2*h) 
        }
        else
        {
            Tmm[j,1] <- 4/h
        }

        
        Tim[j,4] <- 0.25*(t[j]-zz[l-3])*Tmm[j,4] + 0.25*hh2*Tmm[j,3] + 0.25*h3m*Tmm[j,2] + 0.25*h4*Tmm[j,1]
        Tim[j,3] <- 0.25*hht*Tmm[j,3] + h3m*Tmm[j,2]*0.25 + h4*Tmm[j,1]*0.25
        Tim[j,2] <- htm*Tmm[j,2]*0.25 + h4*Tmm[j,1]*0.25
        Tim[j,1] <- ht*Tmm[j,1]*0.25
    }

    return(list(Tmm=Tmm, Tim=Tim))        
}
