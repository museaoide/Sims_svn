normAlmd <- function(Aml, lmdml, A, lmd) {
    nsig <- dim(lmd)[2]
    nv <- dim(lmd)[1]
    Alml <- array(0, c(nv, nv, nsig))
    Al <- Alml
    for (il in 1:nsig) {
        Alml[ , , il] <- exp(-.5 * lmdml[ , il]) * Aml
        Al[ , , il] <- exp(-.5 * lmd[ , il]) * A
    }
    Alml <- matrix(Alml, nv)
    Al <- matrix(Al, nv)
    xp <- Al %*% t(Alml)
    ## algorithm is not idempotent.  Fix.
    ordrng <- vector("numeric", nv)
    crit <- vector("numeric", nv)
    ##---------- 3rd attempt--------------
    ## Make any switch with 1 that increases trace, then any with 2, etc.
    for (iv in 1:nv) {
        for (iv2 in iv:nv) {
            crit[iv2] <- xp[iv2,iv] - xp[iv,iv] + xp[iv,iv2] - xp[iv2,iv2]
        }
        idtr <- which.max(crit[iv:nv]) + iv - 1
        Al[iv:nv, ] <- Al[c(idtr, (iv:nv)[-idtr]), ]
        ordrng[iv] <- idtr
    }
    sf <- diag(Al)
    Al <- (1/sf) * Al
    lmd <- sf * lmd
    return(list(Anormed=Al , lmdnormed=lmd, ordrng=ordrng))
}
