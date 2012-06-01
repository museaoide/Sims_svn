restrictVAR <- function(vout, rmat=NULL, yzrone=NULL, xzrone=NULL, const=NULL) {
  ## restrictions can be specified as rows of rmat, with coefficients applied to elements of By and Bx
  ## stacked as they are in xxi (and then repeated across the equation index), or they can be specified
  ## in yzrone, xzrone.  Each zero element of yzrone or xzrone generates a restriction that sets the corresponding
  ## coefficient in By or Bx to zero.  Both kinds of restrictions can be non-trivial in the same call.
  ncf <- dim(vout$By)[2] * dim(vout$By)[3] + dim(vout$Bx)[2]
  neq <- dim(vout$By)[1]
  ny <- dim(vout$By)[2]
  lags <- dim(vout$By)[3]
  cfarray <- 
  nx <- dim(vout$Bx)[2]
  if (is.null(rmat)) {
    rmat <- matrix(0, 0, ncf *neq)
  }
  if (!is.null(yzrone)) {
    byz <- which(yzrone == 0, arr.ind=TRUE)
    for (ir in 1:dim(byz)[1] ) {
      newrow <- rep(0, neq * ncf)
      newrow[(byz[ir,1] - 1) * ncf + (byz[ir, 3] -1) * ny + byz[ir, 2]] <- 1
      rmat <- rbind(rmat,newrow)
    }
  }
  if (!is.null(xzrone)) {
    bxz <- which(xzrone == 0, arr.ind=TRUE )
    for (ir in 1:dim(bxz) ) {
      newrow <- rep(0,ncf * neq)
      newrow[(byx[ir,1] - 1) * ncf + ny * lags + byx[ir, 2]] <- 1
      rmat <- rbind(rmat, newrow)
    }
  }
  sig <- cov(vout$u)
  rvcv <- rmat %*% kronecker(sig, vout$xxi) %*% t(rmat)
  if(is.null(const)) const <- rep(0, dim(rmat)[1])
  stackedcf <- c(t(cbind(matrix(vout$By, nrow=neq), vout$Bx)))
  gap <- rmat %*% stackedcf - const
  chstat <- t(gap) %*% solve(rvcv) %*% gap
  return(chstat)
}
  
