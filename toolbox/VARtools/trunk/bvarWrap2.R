bvarWrap2 <- function(x, verbose=FALSE) {
    ## For returning detailed results, set verbose=TRUE
    ## Tsigbrk <- invTime(c(1979.75, 1983.0, 2008.0, 2010.0),  slimdata4)
    Tsigbrk <- invTime(c(1979.75, 1983.0, 2008.0),  slimdata4)
    ## here, Tsigbrk is when new sig starts; below we shift it back to be last obs with old sig.
    Lags <- 6
    nv <- 6
    ## enddata <- 2007.75
    enddata <- end(slimdata4)
    T <- dim(window(slimdata4, end=enddata))[1]
    Tsigbrk <- c(0, Tsigbrk - 1)
    Tsigbrk <- c(Tsigbrk[Tsigbrk < T], T)
    nsig <- length(Tsigbrk) - 1
    A <- matrix(0, nv, nv)
    dgx <- seq(1, nv^2, by=nv+1)
    A[dgx] <-  1
    A[-dgx] <- x[1:(nv^2 - nv)]
    lmd <- matrix(0, nv, nsig)
    lmd[ , ] <- x[nv^2 - nv + 1:(nv*nsig)]
    sigfac <- array(0, c(nv, nv, nsig))
    for (isig in 1:nsig) 
        sigfac[ , , isig] <- exp(-.5 * lmd[ , isig]) * t(A)
    sig0 <- crossprod(sigfac[ , , 1])
    vnames <- dimnames(slimdata4)[[2]]
    dimnames(sig0) <- list(vnames, vnames)
    ## --------- set up prior parameters ---------------
    mnprior <- list(tight=5, decay=.5)
    ## vprior <- list(sig=sqrt(diag(sig0)), w=0)
    vprior <- list(sig=c(rep(.01,5),.1), w=0)
    names(vprior$sig) <- dimnames(slimdata4)[[2]]
    urprior <- list(lambda=5, mu=1)
    ybar <- apply(slimdata4[1:6, ], 2, mean, na.rm = TRUE)
    sigfix <- diag(vprior$sig^2)
    ## nstat <- rep(TRUE, 6)
    nstat <- c(rep(TRUE, 4), FALSE, FALSE)
    dimnames(sigfix) <- list(vnames,vnames)
    prior <- varpriorN(nv, nx = 1, lags = Lags, mnprior = mnprior, vprior = vprior, urprior = urprior, ybar = ybar, sigfix=sigfix, nstat=nstat)
    vout <- rfvarKFx(ydata = window(slimdata4, end=enddata), lags = Lags, sigfac = sigfac, Tsigbrk=Tsigbrk, prior = prior)
    lh <- -sum(vout$lh)                 #Note sign flip, for minimization
    attr(lh,"prior") <- list(mnprior=mnprior, vprior=vprior, urprior=urprior, ybar=ybar, sigfix=sigfix, nstat=nstat)
    attr(lh,"sigfac") <- sigfac
    attr(lh, "T") <- T
    attr(lh, "data") <- "slimdata4"
    ## prior on lambda's, to stay away from zeros.
    ## lmscale <- .002
    ## nlmd <- length(c(lmd))
    ## lplmd <- nlmd * (log(2) -  2 * log(lmscale)) + 3 * sum(-lmd) - 3 * sum(log(1 + lmscale^(-2) * exp(-2 * lmd)))
    ## simpler exponential prior on lambdas
    lplmd <- -sum(lmd) - sum(exp(-lmd))  #exp(-lmd) ~ exponential
    ##---------- ev penalty messes up MCMC logic, so omit unless really needed. ---------
    ## ## penalize highly unstable roots
    ## ev <- eigen(sysmat(vout$By))$values
    ## ev <- sum(abs(ev[abs(ev) > 1] - 1))^2*1e3 #(so one root of 1.03 penalized by .9 (weak)
    ##-----------------------------
    ev <- 0
    ##-------------------------------------------------------------------------
    lh <- lh + ev - lplmd                     # correct marginal posterior pdf | lmd, A
    ##-------------------------------------------------------------------------
    ## dsig <- diag(chol(vout$Vb, pivot=TRUE))
    ## if (any(dsig <= 0)) {
    ##     dsig <- -1e20
    ## } else {
    ##     dsig <- sum(log(dsig))
    ## }
    ## lh <- lh + ev - lplmd - dsig        #fixed "-.5 * dsig", 13.7.12
    ## ## dsig component here converts max'd llh into log of mgnl lh given lmd, A.
    attr(lh, "penalty") <- ev
    if(verbose) {
        ## form in-sample residuals---------
        u <- fcastMany(slimdata4, vout$By, vout$Bx, horiz=1:8)
        u1std <- matrix(0, dim(u$u)[1]+Lags, dim(u$u)[3])
        for (isig in 1:nsig) {
            u1std[(Tsigbrk[isig] + 1):Tsigbrk[isig + 1], ] <-
                u$u[(Tsigbrk[isig] + 1):(Tsigbrk[isig + 1]) , 1, ] %*%
                    diag(1/sqrt(diag(crossprod(sigfac[ , , isig]))))
            u1std <- ts(u1std, start=tsp(u$u)[1], freq=tsp(u$u)[3])
        }
        return(list(lh=lh, vout=vout, A=A, lmd=exp(-.5*lmd), llmd = lmd, f=u$fc, u=u$u, u1std=u1std)) #llmd included because exp(-.5*lmd) might lose precision.
    } else {
        return(lh)
    }
}
