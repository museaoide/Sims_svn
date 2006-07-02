HClhObj <- function(parvec, y, x, pyGs, phat0, lags, parPrior, ...) {
  ## To keep this function generic, we do not deduce lags from parvec, but make it a separate
  ## argument.  Only the pij component of parvec is parsed out in HClhObj.
  transMat <- makeTmat(parvec$pij)
  ##  rho <- parvec$rho
  ##  gam <- parvec$gam
  ##  sigsq <- parvec$sigsq
  ##  nState <- length(sigsq)
  nT <- dim(ydata)[1]
  alpha <- parPrior$alpha
  lpijNorm <- sum(lgamma(alpha)) - sum(lgamma(rep(1,nState) %*% alpha)) # (nState columns of pij's)
  if ( is.null(dim(y)) ) {
    dim(y) <- c(length(y),1)
    if (!is.null(x) && dim(x)[1] != dim(y)[1]) stop("x and y must have the same nrow")
  }
  nT <- dim(y)[1]
  psGYfilt <- phat0
  llh <- 0.0
  for (it in (lags+1):nT) {
    lhe <- HClhElement(y=y[it, ], pyGs=pyGs, phat=psGYfilt[it-1, ], transMat=transMat, yl=y[(it-lags):(it-1), ], x=x[it, ], parvec=parvec, )
    psGYfilt <- lhe$phatnew
    llh <- llh + lhe$lhElement
  }
  dObsOut <- setDumObs(parPrior, sigsq, y0, x0) 
  return(llh-lpijNorm)
}
