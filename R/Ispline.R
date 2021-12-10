#' @export
Ispline <- function(y,z)
{
    nz <- length(z)
    z <- sort(z)
    if(any(na.omit(y)<z[1]) | any(na.omit(y)>z[nz])) stop("y en dehors de z")

    ## repeter les noeuds externes
    zz <- rep(0,nz+4)
    zz[1:2] <- rep(z[1],2)
    zz[2+1:nz] <- z
    zz[2+nz+1:2] <- z[nz]

    ## declaration du resultat
    im <- matrix(0,nrow=length(y),ncol=3)
    mm <- matrix(0,nrow=length(y),ncol=3)

    ## boucle sur les y
    for(j in 1:length(y))
    {
        ## encadrer y par les noeuds
        yz <- c(z[-c(1,nz)],y[j])
        iyz <- c(rep(0,nz-2),1)
        l <- 2+which(iyz[order(yz)]==1)

        ## calcul des im et mm

        ht2 <- zz[l+1]-y[j]
        htm <- y[j]-zz[l-1]
        ht <- y[j]-zz[l]
        ht3 <- zz[l+2]-y[j]
        hht <- y[j]-zz[l-2]
        h <- zz[l+1]-zz[l]
        hh <- zz[l+1]-zz[l-1]
        hn <- zz[l+1]-zz[l-2]
        h2n <-zz[l+2]-zz[l-1]
        h2 <- zz[l+2]-zz[l]
        h3 <- zz[l+3]-zz[l]

        if(y[j]<z[nz])
        {
            mm[j,3] <- (3*ht2*ht2)/(hh*h*hn)
            mm[j,2] <- (3*htm*ht2)/(h2n*hh*h)+(3*ht*ht3)/(h2*h*h2n)
            mm[j,1] <- (3*ht*ht)/(h3*h2*h)
        }
        else
        {
            mm[j,1] <- 3/h
        }
                
        im[j,3] <- hht*mm[j,3]/3 + h2n*mm[j,2]/3+h3*mm[j,1]/3
        im[j,2] <- htm*mm[j,2]/3 + h3*mm[j,1]/3
        im[j,1] <- ht*mm[j,1]/3
    }
    
    return(list(im=im,mm=mm))
}
