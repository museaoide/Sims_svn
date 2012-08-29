fcastMany <- function(ydata, By, Bx, xdata=NULL, const=TRUE, horiz) {
  ## calculates forecasts and forecast errors at horizons given by the horiz vector
  ## at each date from lags+1 through end(ydata)-min(horiz).
  ## xdata, if present, must have first dimension exceeding that of ydata by max(horiz)
  ##-----------------------------------------------------------
  lags <- dim(By)[3]
  T <- dim(ydata)[1]
  nv <- dim(By)[1]
  nx <- ifelse( is.null(xdata), 0, dim(xdata)[2])
  nx <- ifelse(const, nx + 1, nx)
  hmax <- horiz[length(horiz)]
  hmin <- horiz[1]
  fc <- array(0, c(dim(ydata)[1] - lags, length(horiz), nv))
  u <- fc
  if (const) xdata <- cbind(xdata, matrix(1, T + hmax, 1))
  for (it in (lags + 1):(T - hmin)) {
    y0 <- ydata[(it - lags):(it-1), ]
    hz <- horiz[(it + horiz) <= T]
    hmax <- max(hz)
    x0 <- xdata[(it - lags):(it - 1 + hmax), ]
    fc[it - lags, 1:length(hz), ] <- fcast(y0, By, Bx, x0, const=FALSE, hmax)[ lags + hz, ]
    ## const=FALSE for fcast because we have already filled x0 with ones.
    ## note that fcast returns initial conditions and forecast all stacked up.
    u[it - lags, 1:length(hz) , ] <- ydata[it + hz, ] - fc[it - lags, 1:length(hz), ]
  }
  dimnames(fc) <- list(NULL, horiz, dimnames(ydata)[[2]])
  dimnames(u) <- dimnames(fc)
  return(list(fc=fc, u=u))
}
