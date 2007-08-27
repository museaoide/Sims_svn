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
  xdum <- t(Matrix(t(xdum), dim(xdum)[2], nv * dim(xdum)[1])) # stacks up vertically nv copies of xdum
  ## Get actual data in regression form
  y <- Matrix(ydata[-(1:lags), , ], ncol=1)
  x <- array(0, c(nT-lags, nv, lags, nc))
  for (il in 1:lags)
    x[ , , il, ] <- ydata[(lags+1):nT - il, , ]
  X <- NULL
  for (ic in 1:nc)
    X <- rbind(X, kronecker(Diagonal(nv), Matrix(x[ , , , ic], nrow = nT - lags))) #each variable in a country gets the same
                                        #matrix of rhs lagged y's
  ## now add all the country-specific constants
  cdum <- array(diag(nv), c(nv, nv, nT-lags))
  cdum <- aperm(cdum, c(3,1,2))
  cdum <- Matrix(cdum, ncol=nv)
  X <- cbind(X, kronecker(Diagonal(nc), cdum))
  nActualRows <- dim(X)[1]
  X <- rbind(X, Matrix(0, nrow=dim(xdum)[1], ncol=dim(X)[2]))
  X[-(1:nActualRows), 1:dim(xdum)[2]] <- xdum
  ##
  ## Now add rows that connect constants to y0's
  y0 <- Matrix(ydata[1:lags, , ], ncol=1) - Matrix(Mu, lags * nv * nc, ncol=1)  #subtracts out the mu from y0=Lambda*cs+mu
  x0 <- Matrix(0, lags*nv*nc, dim(X)[2])
  x0[ , lags * nv + 1:(nc * nv)] <- kronecker(Diagonal(nc), Lambda)
  X <- rbind(X, x0)
  y <- rbind(y,y0)
###!!! still need cs prior
  ## Set up residual covariance matrix
  BigFac <- Matrix(0, length(y), length(y))
  for (ic in 1:nc) {
    icndx <- (ic - 1) * nv * (nT-lags) + 1:(nv * (nT-lags))
    BigFac[icndx, icndx] <- kronecker(SigFac[ , , ic], Diagonal(nT-lags))
  }
  OmFacNdx <- (nT - lags) * nv * ic + 1:(nv * lags * nc)
  stopifnot(max(OmFacNdx) == dim(BigFac)[1]) # delete this line once things are working cleanly
  BigFac[OmFacNdx, OmFacNdx] <- kronecker(Diagonal(nc), OmFac)
  qrX <- qr(BigFac %*% X)
  stdy <- BigFac %*% y
  cf <- qr.coef(qrX, stdy)
  stdresid <- qr.resid(qrX, stdy)
  resid <- y - X %*% cf
  ActualNdx <- c(nT-lags, nv, ic)
  
}
