panelPrior <- function (yinit, xinit, lags, lambda, sig) {
  ## yinit:     Tinit x nvar x ngroup array of initial conditions for y  (Usually Tinit=lags)
  ## xinit:     Tinit x nx x ngroup array of initial country x's.  (maybe just country dummies)
  ## lambda:    overall scale on the co-persistence dummies
  ## sig:       like vprior.sig in rfvar3:  approximate prior modes for std deviations of disturbances
  nvar <- dim(yinit)[2]
  nx <- dim(xinit)[2]
  ngroup <- dim(yinit)[3]
  ybar <- apply(yinit, c(2,3), mean)
  xbar <- apply(xinit, c(2,3), mean)
  ydum <- array(0, c(lags+1, nvar, ngroup))
  xdum <- array(0, c(lags+1, nx, ngroup))
  for (il in 1:(lags+1)) {
    ydum[il, , ] <- ybar
    xdum[il, , ] <- xbar
  }
  for (iv in 1:nvar) {
    ydum[ , iv, ] <- lambda * ydum[ , iv, ] / sig[iv]
  }
  pbreak <- (lags+1) * 1:(ngroup-1)
  ydum <- aperm(ydum, c(1,3,2)) # lags, group, var
  ydum <- matrix(ydum, ncol=nvar)
  return(list(ydum=ydum,pbreak=pbreak))
}
