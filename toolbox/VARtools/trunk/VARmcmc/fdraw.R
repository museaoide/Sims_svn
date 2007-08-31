fdraw <- function(y0, By, Bx, xdata=NULL, const=TRUE, horiz, xxi, smat, df, ndraw) {
  scaledraw <- rvecwish(df, cs=t(smat), ndraw=ndraw)
  nv <- dim(By)[1]
  shocks <- rnorm(nv*nv*horiz*ndraw)
  shocks <- array(shocks, c(nv,nv,horiz,ndraw))
  for (draw in 1:ndraw) {
    shocks <- tensor(scaledraw[ , , draw], shocks[ , , ,draw], 1, 1)
  }
  lags <- dim(y0)[[1]]
  yf <- array(0, nv, horiz, ndraw)  * here need to duplicate logic of fcast, but with shocks, and vectorized over draws
  for (it in 1:horiz)
    yf[ , 
