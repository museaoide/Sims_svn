bvarWrap2 <- function(x, prior, verbose=FALSE) {
    ## For returning detailed results, after convergence, use second return
    ## below and comment out the simple return(lh).
    Tsigbrk <- invTime(c(1979.75, 1983.0, 2008.0, 2010.0),  slimdata4)
    ## here, Tsigbrk is when new sig starts; below we shift it back to be last obs with old sig.
    Lags <- 6
    nv <- 6
    enddata <- 2007.75
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
    vout <- rfvarKFx(ydata = window(slimdata4, end=enddata), lags = Lags, sigfac = sigfac, Tsigbrk=Tsigbrk, prior = prior)
    lh=-sum(vout$lh)
    attr(lh,"prior") <- prior
    attr(lh,"sigfac") <- sigfac
    attr(lh, "T") <- T
    attr(lh, "data") <- "slimdata4"
    ## prior on lambda's, to stay away from zeros.
    lmscale <- .002
    nlmd <- length(c(lmd))
    lplmd <- nlmd * (log(2) -  2 * log(lmscale)) + 3 * sum(lmd) - 3 * sum(log(1 + lmscale^(-2) * exp(2 * lmd)))
    ## penalize highly unstable roots
    ev <- eigen(sysmat(vout$By))$values
    ev <- sum(abs(ev[abs(ev) > 1] - 1))^2*1e3 #(so one root of 1.03 penalized by .9 (weak)
    lh <- lh + ev - lplmd
    attr(lh, "penalty") <- ev
    if(verbose)
        return(list(lh=lh, vout=vout, A=A, lmd=exp(lmd))) # uncomment for analyzing output
    else
        return(lh)
}
