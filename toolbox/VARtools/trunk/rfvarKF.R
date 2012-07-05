rfvarKF <- function(ydata=NA,lags=6,xdata=NULL,const=TRUE,breaks=NULL,lambda=5,mu=2,ic=NULL, sig) {
  ## like rfvar3, but uses KF instead of direct matrix calculations, thereby allows for time-varying
  ## variances in the sig array.
  
