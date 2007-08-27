mleBC <- function(ydata, lags, SigFacs, mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1), mu=10, Lambda, Mu) {
  require(Matrix)
  nT <- dim(ydata)[[1]]
  nv <- dim(ydata)[[2]]
  nc <- dim(ydata)[[3]]
  vp <- varprior(nv=nv,nx=0,lags=lags, mnprior=mnprior,vprior=prior)
  ndobs <- dim(vp$ydum)[1] - (length(vp$pbreaks)+1) * lags
  ydum <- rbind(vp$ydum[vp$breaks, ], vp$ydum[dum(vp$ydum)[1], ])
  ## line above works because each dummy observation takes up exactly lags+1 positions in vp$ydum.  With training
  ## sample component (as in rfvar3 or mgnldnsty), would need more logic.
  xdum <- matrix(vp$ydum[-c(breaks, dim(vp$ydum)[1]), ], nrow=ndobs)
  ## now sum of coeffs dummy obs
  ## Unlike rfvar3.R, here we don't use initial conditions for these.  This might also be better in rfvar3.
  ## Now the observations are scaled by vprior$sig, which means that the mu weight on them must be much bigger here
  ## than in rfvar3.R to have the same strength.  mu of 5 or 10 would be sensible.
  if (mu > 0) {
    ydumsc <- mu * Diagonal(vprior$sig, nrow=nv)
    xdumsc <- array(diag(nv), c(nv, nv, lags))
    xdumsc <- Matrix(aperm(xdumsc, c(1,3,2)), nrow=nv)
    ydum <- rbind(Matrix(ydum), ydumsc)
    xdum <- rbind(Matrix(xdum), xdumsc)
  }
    ## string out these dummy obs, because prior connection of cs to y0 below breaks
    ## I cross X structure of rhs.
  dim(ydum) <- c(dim(ydum)[1]*nv,1)
  xdum <- t(Matrix(t(xdum), dim(xdum)[2], nv * dim(xdum)[1])  # stacks up vertically nv copies of xdum
  
}
