fcast <- function(y0, By, Bx, xdata=NULL, const=TRUE, horiz) {
### By: equations x variables x lags
### Bx: equations x nx
### xdata: (lags+horiz) x nx
### y0: lags x nvar.  Must be a matrix, even if univariate.
### produces horiz-step-ahead forecast.
### Value is:
### yhat: (horiz+lags) x nvar
  if (is.null(dim(y0)))
    lags <- length(y0)
  else
    lags <- dim(y0)[1]
  stopifnot( lags == dim(By)[3] )
  stopifnot( is.null(xdata) || (is.null(dim(xdata)) && length(xdata) == horiz+lags) || horiz+lags == dim(xdata)[1] )
  if (const) {
    if (is.null(xdata))
      xdata <- matrix(1,horiz+lags,1)
    else
      xdata <- cbind(xdata, matrix(1, horiz+lags, 1))
  }
  if (is.null(dim(y0))) dim(y0) <- c(length(y0),1)
  nvar <- dim(y0)[2]
  if (is.null(dim(xdata))) xdata <- matrix(xdata,ncol=1)
  if (is.null(dim(Bx))) Bx <- matrix(Bx,ncol=1)
  nx <- dim(Bx)[2]
  yhat <- matrix(0,horiz+lags,nvar)
  yhat[1:lags,] <- y0
  Bmat <- aperm(By,c(3,2,1))
  Bmat <- Bmat[seq(lags,1,by=-1),,]
  dim(Bmat) <- c(lags*nvar,nvar)
  for (it in 1:horiz){
    ydata <- yhat[it:(it+lags-1),,drop=FALSE]
    yhat[lags+it,] <- apply(Bmat*matrix(ydata,dim(Bmat)[1],dim(Bmat)[2]),2,sum)+xdata[lags+it,]%*%t(Bx)
  }
  if(!is.null(dimnames(y0))){
    dimnames(yhat) <- list(NULL,dimnames(y0)[[2]])
  }
  if(is.ts(y0)){
    yhat <- ts(yhat,start=tsp(y0)[1],freq=tsp(y0)[3])
  }
  return(yhat)
}
