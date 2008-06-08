SVARlh0 <- function(pvec, idmat, sigma, T) {
  ## idmat is a logical matrix, TRUE where A0 has a non-zero coefficient
  ## pvec  is the vector of parameter values that fill A0[idmat].
  ## sigma is the covariance matrix of estimated residuals from the reduced
  ##       form VAR
  ## T     is the sample size.
  ## This function returns minus the likelihood, so it can be used directly in csminwel
  n <- dim(idmat)[1] # better be same as dim(idmat[2])
  A0 <- matrix(0,n,n)
  A0[idmat] <- pvec
  lh <- SVARlh(A0, sigma, T)
  ## penalty for contemporaneous response of FF (variable 1) to M (variable 6) negative
  ## and contemp policy response to pcrm small
  Tdum <- 20
  lh <- lh - Tdum*log(2*pi) + 2*Tdum*log(A0[1,1]) - Tdum*((.005*A0[1,1] + .007*A0[1,6])^2/2 +(.005*A0[1,1] + .05 * A0[1,7])^2/2)
  return(-lh)
}
