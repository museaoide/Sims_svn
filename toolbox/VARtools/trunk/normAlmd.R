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
    ## Algorithm tries reordering up to nv times to find an invariant ordering,
    ## then gives up and returns nv'th reordering and noloop=FALSE
    ordrng <- 1:nv
    crit <- vector("numeric", nv)
    Alold <- Al
    noloop <- 0
    for (ntrial in 1:nv) {
        ## Make any switch with 1 that increases trace, then any with 2, etc.
        for (iv in 1:nv) {
            for (iv2 in iv:nv) {
                crit[iv2] <- xp[iv2,iv] - xp[iv,iv] + xp[iv,iv2] - xp[iv2,iv2]
            }
            idtr <- which.max(crit[iv:nv])
            ordrng[iv:nv] <- ordrng[iv:nv][c(idtr , (1:(nv-iv+1))[-idtr])]
            Al <- Al[ordrng, ]
            xp <- xp[ordrng, ]
        }
        if (isTRUE(all.equal(Al, Alold))) {
            noloop <- ntrial
            break
            Alold <- Al
        }
    }
    sf <- diag(Al)
    Al <- (1/sf) * Al
    lmd <- sf * lmd
    return(list(Anormed=Al , lmdnormed=lmd, ordrng=ordrng, noloop=noloop))
}
