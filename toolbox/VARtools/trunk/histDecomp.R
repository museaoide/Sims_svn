histDecomp <- function(vout, vts, xdata=NULL, const=TRUE, orthmat=NULL) {
  if(is.null(orthmat)) orthmat <- chol(cov(vout$u))
  lags <- dim(vout$By)[3]
  nv <- dim(vout$By)[1]
  T <- dim(vts)[1]
  horiz <- T - lags
  if (is.null(xdata)) nx <- 0 else nx <- dim(xdata)[2]
  if (const) nx <- nx + 1
  ## form deterministic part
  ydet <- fcast(vts[1:lags, ], vout$By, vout$Bx, xdata, const, horiz)
  ydecomp <- array(0, c(nv, nv + 1, T))
  ydecomp[ , nv+1, ] <- t(ydet)
  orthmati <- solve(orthmat)
  orthshock <- vout$u[(lags + 1):T, ] %*% orthmati
  for (iv in 1 : nv) {
    shocks <- matrix(0, T - lags, nv)
    shocks[ , iv] <- orthshock[ , iv]
    shocks <- shocks %*% orthmat
    ydecomp[ , iv, ] <-  t(fcast(matrix(0, lags, nv), vout$By, vout$Bx, xdata=NULL, const=FALSE, horiz, shocks))
  }
  ydecompStacked <- apply(ydecomp, c(1,3), cumsum)
  ydecompStacked <- aperm(ydecompStacked, c(2, 1, 3))
  dn2 <- dimnames(vout$By)[[1]]
  if (const) dn2 <- c(dn2, "const")
  dimnames(ydecomp) <- list(dimnames(vout$By)[[1]], dn2, NULL)
  dimnames(ydecompStacked) <- dimnames(ydecomp)
  return(list(ydec=ydecomp, ydecStack=ydecompStacked))
}
  
