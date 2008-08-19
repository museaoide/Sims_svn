SigInit <- function(A, Omega, T, mu0, Sig0, Tfac=1, ct=FALSE) {
  ## The system is y(t) = A(t) %*% y(t-1) + eps(t) in discrete time, ydot = A %*% y + eps in
  ## continuous time, with var(eps(t))=Omega.
  ## If the system has a constant term, it is assumed to be represented
  ## by a row of A that has in discrete time only zeros and a single 1 entry, and a corresponding
  ## row and column of Omega that are all zero.  In continuous time the row of A is all zeros.
  ## T:     sample size.
  ## mu0:   usually a vector of 0's, except for a 1 corresponding to the constant.
  ## Sig0:  A big covariance matrix, to which the distribution of the initial conditions
  ##        converges as roots approach the boundary of the stationary region.
  ##        It should imply standard errors several times as large as the y data itself.
  ##        It should be zero for any constant or deterministic trend components of y.
  ##        If big enough, it should have little effect on the posterior density, but
  ##        it can in principle have important effects on model comparisons, so should be
  ##        varied to check sensitivity in model comparisons.  
  ## Tfac:  Tfac*T is the minimum approximate half-life of a stationary component that
  ##        will be treated as possibly non-stationary, with the weight on non-stationarity
  ##        increasing as the actual half-life increases.
  ## ct:    Is this a continuous time model (so Re() rather than abs() ranks roots)?
  ##
  if (!is.loaded("zhseqr")) dyn.load("usr/lib/liblapack.so")
  wtfcn <- function(x) max(1 - x + sin(2*pi*x)/(2*pi),0)
  if (ct) {
    div <- c(-1/(Tfac*T), -1000*.Machine$double.eps) # might need to adjust div[2]
  } else {
    div <- c(1-1/(Tfac*T), 1-1000*.Machine$double.eps)
  }
  sca <- blkDglz(A, div, ctOrder=ct)
  nlow <- sca$blockDims[1]
  if (nlow == 0) {
    mf1 <- rep(0,0)
    v1 <- matrix(1,0,0)
    lowx <- NULL
  } else {
    lowx <-  1:nlow 
    ## mf1 <- solve(diag(nlow) - sca$D[lowx,lowx], sca$Pinv[lowx, ] %*% mu0)
    mf1 <- matrix(0,nlow,1)
    omega1 <- sca$Pinv[lowx, ] %*% Omega %*% t(Conj(sca$Pinv[lowx, , drop=FALSE]))
    if (ct) {
      v1 <- sylvester(-sca$D[lowx, lowx], omega1)$X
    } else {   
      v1 <- doubling(sca$D[1:nlow,1:nlow], omega1)
    }
  } 
  nmid <- sca$blockDims[2]
  if (nmid == 0) {
    v2 <- matrix(1,0,0)
    mf2 <- rep(1,0)
    midx <- NULL
  } else {
    midx <- (nlow+1):(nlow + nmid)
    mf2ns <- sca$Pinv[midx, ]  %*% mu0  # mean as if non-stationary
    ## mf2s <- solve( diag(nmid) - sca$D[midx, midx], sca$Pinv[midx, ] %*% mu0) #mean as if stationary
    mf2s <- matrix(0,nmid,1)
    if (ct) {
      wt <- wtfcn((Re(diag(sca$D[midx, midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    } else {
      wt <- wtfcn((abs(diag(sca$D[midx,midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    }
    mf2 <- mf2s * wt + mf2ns * (1-wt) #wt is the weight on stationary model.
    omega2 <- sca$Pinv[midx, ] %*% Omega %*% t(Conj(sca$Pinv[midx, , drop=FALSE]))
    if (ct) {
      v2 <- sylvester(-sca$D[midx, midx], omega2)$X
    } else {   
      v2 <- doubling(sca$D[midx,midx], omega2)
    }
    v20 <- sca$Pinv[midx, ] %*% Sig0 %*% t(Conj(sca$Pinv[midx, , drop=FALSE]))
    wta <- c(sqrt(wt))               # c() to strip dimension attribute
    wtb <- c(sqrt(max(1-wt,0)))
    v2 <- (wta %o% wta) * v2 + (wtb %o% wtb) * v20
  }
  nhi <- sca$blockDims[3]
  if (nhi == 0) {
    mf3 <- rep(0,0)
    v3 <- matrix(1,0,0)
    hix <- NULL
  } else {
    hix <- (nlow + nmid + 1):(nlow + nmid + nhi)
    mf3 <- sca$Pinv[hix, ] %*% mu0
    v3 <- sca$Pinv[hix, ] %*% Sig0 %*% t(Conj(sca$Pinv[hix, , drop=FALSE]))
  }
  muout <- sca$P %*% c(mf1, mf2, mf3)
  vout <- matrix(0, nlow + nmid + nhi, nlow + nmid + nhi)
  vout[lowx, lowx] <- v1
  vout[midx, midx] <- v2
  vout[hix, hix] <- v3
  vout <- sca$P %*% vout %*% t(Conj(sca$P))
  return(list(mu=muout , v=vout, P=sca$P, Pinv=sca$Pinv))
}
