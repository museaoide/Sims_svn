normAlmd <- function(Aml, A, lmd) {
    ## we assume columns of Aml have already been scaled so all(diag(crossprod(Aml)) == 1).
    ## For this to be reliable, variances of r.f. innovations must be same
    ## order of magnitude across variables.
    
    sfa <- sqrt(apply(A^2, 2, sum))
    A <- A %*% diag(1/sfa)
    cvmat <- abs(crossprod(Aml, A))
    neq <- dim(A)[2]
    ordrng <- vector("numeric", neq)
    ordrng[1] <- which.max(cvmat[1, ])
    for (iq in 2:neq) {
        ndx <- which.max(cvmat[iq, -ordrng[1:(iq-1)]])
        ordrng[iq] <- setdiff(1:neq, ordrng[1:(iq-1)])[ndx]
    }
    A <- A[ , ordrng]
    sfa <- diag(A)
    A <- A %*% diag(1/sfa)          # So i'th equation has unit coeff on i'th vble.
    lmd <- sfa * lmd[ordrng, ]    # note:  equiv to sfa^2 %*% lmd[ordrng, ]
    return(list(Anormed=A , lmdnormed=lmd, ordrng))
}
