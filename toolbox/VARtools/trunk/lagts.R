lagts <- function(xts,lags) {
  ## The general case is that xts is an mts object and lags a list of length equal to the
  ## number of series in xts.  The lags specified in lags are added to the xts object as 
  ## additional columns, and the resulting mts object is trimmed to eliminate NA's that would 
  ## otherwise arise. The series in the returne object retain their original names, with integers 
  ## appended indicating how far they have been lagged.
  ## In experimenting with regression models with vayring lag lengths, it is a good idea to set up
  ## a single mts object with this program at the start, so that the same sample size is being
  ## used in comparing models with different lags.
  ## Note that the series are trimmed to eliminate all NA's, so that an NA in the middle of any
  ## series in xts  will result in extreme truncation.
  ## xts can also be a single time series.  lags can be just a single vector of lags, in which case
  ## all series in xts are included with the same lags.
  ##---------------------
  xtsout <- xts
  vnames <- if(is.null(dim(xts))) "y" else dimnames(xts)[[2]]
  nv <- if(!is.null(dim(xts))) dim(xts)[2] else 1
  if (!is.list(lags) || length(lags) ==1) {
    lags <- unlist(lags)
    for (il in lags) {
      xtsout <- cbind(xtsout, lag(xts, -il))
    }
    ## group lags of the same variable
    ndx <- matrix(1:(nv*lags), nv, lags)
    xtsout <- xtsout[ , c(t(ndx))]
    dimnames(xtsout)[[2]][1:nv] <- vnames
    for (iv in 2:nv) {
      dimnames(xtsout)[[2]][(iv-1)*nv + (1:lags)] <- paste(vnames[iv], 1:lags, sep="")
    }
  } else {
    for (iv in 1:nv) {
      if (!is.null(lags[[iv]])) {
        for (il in lags[[iv]]) {
          xtsout <- cbind(xtsout, lag(xts[ , iv], -il))
        }
      }
    }
    namecount <- nv
    dimnames(xtsout)[[2]][1:nv] <- vnames
    for (iv in 1:nv) {
      if (!is.null(lags[[iv]])) {
        lvec <- lags[[iv]]
        dimnames(xtsout)[[2]][(namecount + 1):(namecount + length(lvec))] <- 
          paste(vnames[iv], lvec, sep="")
        namecount <- namecount + length(lvec)
      }
    }
  }
  xtsout <- trimts(xtsout)
  return(xtsout)
}