mleBC <- function(ydata, lags, SigFacs, mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1), mu=10, Lambda, Mu, cbar, csigFac=NULL) {
  require(Matrix)
  nT <- dim(ydata)[[1]]
  nv <- dim(ydata)[[2]]
  nc <- dim(ydata)[[3]]
  vp <- varprior(nv=nv,nx=0,lags=lags, mnprior=mnprior,vprior=vprior)
  ndobs <- dim(vp$ydum)[1] - (length(vp$pbreaks)+1) * lags
  ydum <- rbind(vp$ydum[vp$pbreaks, ], vp$ydum[dim(vp$ydum)[1], ])
  ## line above works because each dummy observation takes up exactly lags+1 positions in vp$ydum.  With training
  ## sample component (as in rfvar3 or mgnldnsty), would need more logic.
  xdum <- matrix(vp$ydum[-c(vp$pbreaks, dim(vp$ydum)[1]), ], nrow=ndobs)
  ## now sum of coeffs dummy obs
  ## Unlike rfvar3.R, here we don't use initial conditions for these.  This might also be better in rfvar3.
  ## Now the observations are scaled by vprior$sig, which means that the mu weight on them must be much bigger here
  ## than in rfvar3.R to have the same strength.  mu of 5 or 10 would be sensible.
  if (mu > 0) {
    ydumsc <- mu * Diagonal(x=vprior$sig, n=nv)
    xdumsc <- array(diag(nv), c(nv, nv, lags))
    xdumsc <- Matrix(aperm(xdumsc, c(1,3,2)), nrow=nv)
    ydum <- rBind(Matrix(ydum), Matrix(ydumsc))
    xdum <- rBind(Matrix(xdum), xdumsc)
  }
  ## string out these dummy obs, because prior connection of cs to y0 below breaks
  ## I cross X structure of rhs.
  dim(ydum) <- c(dim(ydum)[1]*nv,1)
  xdum <- t(Matrix(as.vector(t(xdum)), dim(xdum)[2], nv * dim(xdum)[1])) # stacks up vertically nv copies of xdum
  ## Get actual data in regression form
  y <- Matrix(ydata[-(1:lags), , ], ncol=1)
  x <- array(0, c(nT-lags, nv, lags, nc))
  for (il in 1:lags)
    x[ , , il, ] <- ydata[(lags+1):nT - il, , ]
  X <- Matrix(0, nrow=0, ncol = nv * nv * lags)
  for (ic in 1:nc)
    X <- rBind(X, kronecker(Diagonal(nv), Matrix(x[ , , , ic], nrow = nT - lags))) #each variable in a country gets the same
                                        #matrix of rhs lagged y's
  ## now add all the country-specific constants
  cdum <- array(diag(nv), c(nv, nv, nT-lags))
  cdum <- aperm(cdum, c(3,1,2))
  cdum <- Matrix(cdum, ncol=nv)
  X <- cBind(X, kronecker(Diagonal(nc), cdum))
  nActualRows <- dim(X)[1]
  X <- rBind(X, Matrix(0, nrow=dim(xdum)[1], ncol=dim(X)[2]))
  X[-(1:nActualRows), 1:dim(xdum)[2]] <- xdum
  ##
  ## Now add rows that connect constants to y0's
  y0 <- Matrix(ydata[1:lags, , ], ncol=1) - Matrix(Mu, lags * nv * nc, ncol=1) #subtracts out the mu from y0=Lambda*cs+mu
  x0 <- Matrix(0, lags*nv*nc, dim(X)[2])
  x0[ , lags * nv *nv + 1:(nc * nv)] <- kronecker(Diagonal(nc), Lambda)
  X <- rBind(X, x0)
  y <- rBind(y,y0)
  ## prior on country-variable constants
  if (!is.null(csigFac)) {
    y <- rBind(y, cbar)
    xc <- Matrix(0, nv*nc, dim(X)[2])
    xc[ , dim(X)[2] - nv * nc + 1:(nv * nc)] <- Diagonal(nv*nc)
  }
  ## Set up residual covariance matrix
  BigFac <- Matrix(0, length(y), length(y))
  for (ic in 1:nc) {
    icndx <- (ic - 1) * nv * (nT-lags) + 1:(nv * (nT-lags))
    BigFac[icndx, icndx] <- kronecker(SigFac[ , , ic], Diagonal(nT-lags))
  }
  OmNdx <- (nT - lags) * nv * ic + 1:(nv * lags * nc)
  BigFac[OmNdx, OmNdx] <- kronecker(Diagonal(nc), OmFac)
  if (!is.null(csigFac)) {
    csigNdx <- max(OmNdx) + 1:(nv * nc)
  }
  BigFac[csigNdx, csigNdx] <- kronecker(Diagonal(nc), csigFac)
  stopifnot (max(csigNdx) == dim(BigFac)[2])  #This line can come out once things are running smoothly
  qrX <- qr(BigFac %*% X)
  stdy <- BigFac %*% y
  cf <- qr.coef(qrX, stdy)
  stdresid <- qr.resid(qrX, stdy)
  resid <- y - X %*% cf
  By <- array(cf[1:(nv * lags * nv)], c(nv, lags, nv)) #variables by lags by equations (order of X columns)
  By <- aperm(By, c(3,1,2))             #now equations by variables by lags, to match other VAR programs
  cs <- Matrix(cf[-(1:(nv * lags * nv))], nv, nc)
  cs <- t(cs)
  dataNdx <- 1:(nv * (nT - lags) * nc)
  vpNdx <-  dataNdx[2] + 1:ndobs
  scNdx <- vpNdx[2] + 1:(nv * nv)
  y0Ndx <- OmNdx
  cspriorNdx <- csigNdx
  return(list(By=By, cs=cs, resid.std.data=stdresid[dataNdx], resid.raw.data=resid[dataNdx], resid.vp=stdresid[vpNdx],
              resid.sc=stdresid[scNdx], resid.std.y0=stdresid[y0Ndx], resid.raw.y0=resid[y0Ndx], resid.cs=stdresid[cspriorNdx]))
}
