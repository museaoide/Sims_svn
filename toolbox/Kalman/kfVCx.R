kfVCx <- function(y, X, shat, sig, M) {
    ## s is the state, and the plant equation is s(t)=G %*% s(t-1)+t(M) %*% e, where e is
    ## N(0,I).  The observation equation is y(t)=Hs(t).  The prior distribution for
    ## s is N(shat,sig).  To handle the standard Kalman Filter setup where the observation
    ## equation is y(t)=Hs(t)+Nu(t), expand s to become [s;v] (where v=Nu), expand H to [H I], replace
    ## G by [G 0;0 0], and replace M with [M 0;0 N].  The resulting model has no "error" in the
    ## observation equation but is equivalent to the original model with error in the state equation.
    ## The posterior on the new state is N(shatnew,signew) and lh is a two-dimensional vector containing
    ## the increments to the two component terms of the log likelihood function.  They are added 
    ## to form the log likelihood, but are used separately in constructing a concentrated or marginal
    ## likelihood. fcsterr is the error in the forecast of y based on the prior.
    ## ----------------
    ## This version specializes to the case of a constant-coefficient VAR, where H= cbind(kronecker(I, X[it, ], I)
    ## and G is the identity (for the constant coefficients) with zeros (for the shocks in the state vector)
    ## appended in the lower right.  Also, kf2's M is all zeros except for the lower right nvar x nvar.  By
    ## bypassing lots of multiplication by zero, this version is faster than generic kf2by a factor of four for
    ## a 7-variable, 13-lag VAR.  Here M is the
    ## cholesky factor of the VAR shock covariances, not the full state equation M (whichis full of
    ## zeros).  With M constant, plain rfvar3 is more efficient.
    ## ----------------------
    ##cat("\nstart kfVCx", proc.time())
    require(tensor)
    SMALLSV <- 1e-30
    if (!is.null(dim(y))) {
        nv <- dim(y)[2]
        nx <- dim(X)[2]
        T <- dim(y)[1]
    } else {
        nv <- length(y)
        nx <- length(x)
        T <- 1
        y <- matrix(T, nv)
        X <- matrix(T, nx)
    }
    nXX <- nv*nx
    nstate <- length(shat)
    omega <- array(sig, c(nx, nv, nx, nv))
    sigObs <- crossprod(M)
    stopifnot(nstate == length(shat))
    ## ho: kron(I,X) %*% sighat
    ## ho <- array(0, c(T, nv, nx, nv))
    ## sv <- array(0, c(T, nv, T, nv))
    ##cat("\nstart double X prod loop", proc.time())
    ## for (iv in 1:nv) {
    ##     ho[ , iv, , ] <- tensor(X, omega[ , iv, ,  ], 2, 1)
    ##     for (jv in 1:nv)
    ##         sv[ , iv, , jv] <- tensor(ho[ , iv, , jv ], X, 2, 2) + sigObs[iv, jv] * diag(T)
    ## }
    ho <- tensor(X, omega, 2, 1)  #T, nv, nx, nv
    sv <- tensor(ho, X, 3, 2)     #T, nv, nv, T
    sv <- aperm(sv, c(1,2,4,3))
    for (irow in 1:nv)
        for (icol in 1:nv)
            sv[ , irow, , icol] <- sv[ , irow, , icol] + sigObs[irow, icol] * diag(T)
    ##svdsv <- svd( matrix(sv, nv * T, nv * T ))
    ##cat("\nstart qr", proc.time)
    qrsv <- qr( matrix(sv, nv * T, nv * T ), LAPACK=TRUE)
    rankr <- match(TRUE, abs(diag(qrsv$qr)) < SMALLSV)
    if (is.na(rankr)) {
        rankr <- nv * T
    } else {
        rankr <- rankr - 1
        qrsv$qr <- qrsv$qr[ , 1:rankr]
        qrsv$aux <- qrsv$aux[1:rankr]
    }
    shatmat <- matrix(shat,nx, nv )
    fcsterr <- y - X %*% shatmat  
    lh <- c(0,0)
    ## we're relying on qr( ,LAPACK=TRUE) always putting zeros in diagonal of R on bottom.
    ## Note that in the case of 0 < rankr < nv * T, we have to account for possible pivot
    ## in QR
    if (rankr == 0) { # Observation is uninformative. Propagate state.
        shatnew <- shat
        signew <- omega
        if (!all(abs(fcsterr) < 1e-7)) warning("Uninformative H but non-zero fcsterr")
    } else if (rankr < nv *T) {
        R <- qr.R[1:rankr, ]
        d <- diag(R)
        R <- R[ , sort.list(qrsv$pivot)]
        Q <- qr.Q(arsv)[ , 1:rankr]
        RQ <- R %*% Q
        oie <- Q %*% solve(RQ, crossprod(Q, c(fcsterr)))
        lh[1] <- -.5 * c(fcsterr) %*% oie
        lh[2] <- -.5 * sum(log(abs(d)))
        ho <- matrix(ho, T*nv, nx*nv)
        shatnew <- shat + oie %*% ho   # works because oie is a vector, not a matrix
        qho <- crossprod(Q, ho)
        signew <- sig - crossprod(qho , solve(RQ, qho)) 
    } else {
        ## first0 <- match(TRUE, svdsv$d < SMALLSV)
        ## if (is.na(first0)) first0 <- nv * T + 1
        ## u <- svdsv$u[ , 1:(first0-1), drop=FALSE]
        ## v <- svdsv$v[ , 1:(first0-1), drop=FALSE]
        ## d <- svdsv$d[1:(first0-1), drop=FALSE]
        ## svifac <- (1/sqrt(d)) * t(u)   #diag(1/sqrt(d)) %*% t(u)
        ## ferr <- svifac %*% c(fcsterr)
        ##lh[1] <- -.5 * crossprod(ferr)
        oie <- solve(qrsv, c(fcsterr))
        lh[1] <- -.5 * c(fcsterr) %*% oie
        lh[2] <- -.5 * sum( log(abs(diag(qr.R(qrsv)))) )
        ## svoifac <-svifac %*% matrix(ho, nv * T, nx * nv)
        ## shatnew <- crossprod(svoifac, ferr) + shat
        ## signew <- sig - crossprod(svoifac)
        ho <- matrix(ho, T*nv, nx*nv)
        shatnew <- shat + oie %*% ho
        signew <- sig - crossprod(ho, solve(qrsv, ho))
    }
    return(list(shat=shatnew, sig=signew, lh=lh, fcsterr=fcsterr))
}
