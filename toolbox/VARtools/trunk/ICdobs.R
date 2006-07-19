ICdobs <- function(y0, x0, horizon, Ay, Ax, filter, ywt, target ) {
  ## dummy observation component of log likelihood, asserting prior beliefs about forecasts from initial
  ## conditions.  If yhat are forecasts from IC's, this function returns the sum of squares of
  ## convolve(filter, yhat-target)*ywt.  E.g. filter is first difference, ywt is t^2*W, where W is a nvar
  ## by nvar matrix, target=0,  which would tend to flatten forecasts, especially at distant horizons.
  ## With the filter a second difference, forecasts just get pushed toward linear trend.
  ##
  ny <- dim(Ay)[1]; nx <- dim(Ax)[2]; lag <- dim(Ay)[3]
  stopifnot(dim(Ay)[2] == ny, dim(Ax)[1] == ny)
  if (!is.null(ywt)) {
    screp <- rep(0,ndo)
    ndo <- dim(ywt)[4]
    stopifnot( identical(dim(Ay),dim(ywt)[1:3]),  ndo==length(rhs),  identical(dim(xwt)[1:2], dim(Ax)))
    for (i in 1:ndo) {
      screp[i] <- sum(Ay * ywt[,,,i]) + sum(Ax * xwt[,,i]) - rhs[i]
    }
  }else {
    screp <- 0
  }
  if (!is.null(ywti)) {
    ndoi <- dim(ywti)[4]
    hrz <- dim(ywti)[3]
    stopifnot( ndoi==length(rhsi) )
    screpi <- rep(0,ndoi)
    Ayi <- impulsdt(Ay,hrz)
    for (i in 1:ndoi) {
      screpi[i] <- sum(Ayi * ywti[,,,i]) - rhsi[i]
    }
  } else {
    screpi <- 0
  }
  return(-.5 * sum(screp^2,screpi^2))
}
