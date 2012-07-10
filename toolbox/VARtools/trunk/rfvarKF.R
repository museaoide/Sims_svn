rfvarKF <- function(ydata=NA,lags=6,xdata=NULL,const=TRUE,breaks=NULL, sigfac, prior) {
  ## ydata:    T x n data matrix, preferably an mts object
  ## xdata:    T x k exogenous data matrix.  Omit if only a constant is needed and const=-TRUE.
  ## breaks:   breaks in the data.  The first lags data points after a break are used
  ##           as new initial conditions, not data points for the fit.  Breaks separated by less than lags
  ##           result in omission of stretches of data.
  ## sigfac:  T x n x n array of sqrts of covariance matrices for disturbances.  crossprod(sigfac) is residual variance
  ## prior:    nkl-dimensional prior mean and covariance matrix for coefficients (as prior$mean and prior$vcv)
    if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
    T <-dim(ydata)[1]
    nvar<-dim(ydata)[2]
    ##nox=isempty(xdata)
    if (const) {
      xdata <- cbind(xdata,matrix(1,T,1))
    }
    nox <- identical(xdata,NULL)
    if(!nox){
      T2 <- dim(xdata)[1]
      nx <- dim(xdata)[2]
    } else {
      T2 <- T; nx <- 0; xdata<- matrix(0,T2,0)
    } 
    ## note that x must be same length as y, even though first part of x will not be used.
    ## This is so that the lags parameter can be changed without reshaping the xdata matrix.
    ## ------------------------
    if (!identical(T2,T)) {
      print('Mismatch of x and y data lengths')
      return()
    }
    if (identical(breaks,NULL))
      nbreaks <- 0
    else {
      nbreaks<-length(breaks)
    }
    breaks <- c(0,breaks,T)
    if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
      stop("list of breaks must be in increasing order\n")
    smpl <- NULL
    for (nb in 2:(nbreaks + 2)) {
      if ( breaks[nb] > breaks[nb-1] + lags )
        smpl <- c(smpl, (breaks[nb-1] + lags + 1):breaks[nb])
    }
    ## With logic above, one can use an mts-type ydata and omit sections of it by including sequences of breaks separated by
    ## less than lags+1.  E.g. with lags=6, monthly data, breaks=rbind(c(1979,8), c(1980,2), c(1980,8), c(1980,12)) omits
    ## Sep 1979 through Dec 1981, plus 6 months after that, which are initial conditions for the next sample segment.
    shat <- prior$mean
    sighat <- prior$variance
    for (ib in 2:(nbreaks+1)) {
      kfout <- kf2(y[(breaks[ib+1]:
