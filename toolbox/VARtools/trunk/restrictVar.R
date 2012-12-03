restrictVAR <- function(vout, type=c("3", "KF"), rmat=NULL, yzrone=NULL, xzrone=NULL,
                        const=NULL, cyzr=NULL, cxzr=NULL) {
  ## restrictions can be specified as rows of rmat, with coefficients applied to elements of By and Bx
  ## stacked as they are in xxi (and then repeated across the equation index), or they can be specified
  ## in yzrone, xzrone.  Each zero element of yzrone or xzrone generates a restriction that sets the corresponding
  ## coefficient in By or Bx to zero (or to a constant, if !is.null(const)).  Both kinds of restrictions
  ## can be non-trivial in the same call.
  #-------------------------------------------
  ## type:     vout as from rfvar3 ("3") or as from rfvarKF ("KF")
  ## const:    the right hand side of rmat %*% coeff = const, not the constant in the var.
  ## cyzr, cxar:  If using yzrone, xzrone with non-trivial constants, leave const=NULL and specify
  ##           constants with cyzr and cxzr
  ##------------------------------------------
  ## sc:       The Schwarz criterion rejects the restriction if the chisq value plus the sc value
  ##           is positive.  This version of the sc is scale-sensitive.  Variables with higher
  ##           variance are penalized more strongly, as with a prior that expects higher-variance
  ##           variables to explain more variance.
  ##
  ##------------------------------------------------
  ncf <- dim(vout$By)[2] * dim(vout$By)[3] + dim(vout$Bx)[2]
  neq <- dim(vout$By)[1]
  ny <- dim(vout$By)[2]
  lags <- dim(vout$By)[3]
  nx <- dim(vout$Bx)[2]
  if (is.null(rmat)) {
    rmat <- matrix(0, 0, ncf *neq)
  }
  if (!is.null(yzrone)) {
    byz <- which(yzrone == 0, arr.ind=TRUE)
    nrstr <- dim(byz)[1]
    if (is.null( cyzr)) cyzr <- array(0, dim(yzrone))
    for (ir in 1:nrstr ) {
      newrow <- rep(0, neq * ncf)
      newrow[(byz[ir,1] - 1) * ncf + (byz[ir, 3] -1) * ny + byz[ir, 2]] <- 1
      rmat <- rbind(rmat,newrow)
    }
    const <- c(const, cyzr[byz])
  }
  if (!is.null(xzrone)) {
    bxz <- which(xzrone == 0, arr.ind=TRUE )
    nrstr <- dim(bxz)[1]
    if (is.null(cxzr)) cxzr <- matrix(0, neq, nx)
    for (ir in 1:nrstr)  {
      newrow <- rep(0,ncf * neq)
      newrow[(bxz[ir,1] - 1) * ncf + ny * lags + bxz[ir, 2]] <- 1
      rmat <- rbind(rmat, newrow)
    }
    const <- c(const, cxzr[bxz])
  }
  svdr <- svd(rmat)
  if (max(abs(svdr$d)) > 1e10 * min(abs(svdr$d))){
    error("restrictions not full rank")
  }
  ## Note that t(rv) spans the same space as rmat, so the restrictiosn are crossprod(v,coeffs)=gamma
  rv <- svdr$v
  if (type == "3") {
    sig <- cov(vout$u)
    svdsig <- svd(sig)
    singsig <- (max(svdsig$d) > 1e10 * min(svdsig$d))
    svdxxi <- svd(vout$xxi)
    singxxi <- (max(svdxxi$d) > 1e10 * min(svdxxi$d))
    singv <- singsig || singxxi
    if(!singv) {
          ## schwarz <- rmat %*% kronecker(svdsig$u %*% diag(1/sqrt(svdsig$d)), svdxxi$u %*% diag(1/sqrt(svdxxi$d)))
          schwarz <- kronecker((1/sqrt(svdsig$d)) * t(svdsig$u), (1/sqrt(svdxxi$d)) * t(svdxxi$u)) %*% rv
        }
  } else {                              #type=="KF"
    svdVb <- svd(vout$Vb)
    ## schwarz <- rmat %*% svdVb$u %*% diag(1/sqrt(svdVb$d)) #below is more efficient version for large Vb
    schwarz <- (1/sqrt(svdVb$d)) * (t(svdVb$u) %*% rv)
  }
  ## T <- if (type == "3") dim(vout$u)[1] else dim(vout$fcsterr)[1]
  df <- dim(rmat)[1]
  schwarz <- -2 * sum(log(diag(chol(crossprod(schwarz)))))   +
     df * log(2 * pi)
  if(is.null(const)) const <- rep(0, dim(rmat)[1])
  stackedcf <- c(t(cbind(matrix(vout$By, nrow=neq), vout$Bx)))
  gap <- rmat %*% stackedcf - const
  svdv <- svd(rmat %*% vout$Vb %*% t(rmat))
  chstat <- (1/sqrt(svdv$d)) * (t(svdv$u) %*% gap)
  chstat <- crossprod(chstat)
  return(list(chiSquared=chstat, df=df, sc=schwarz))
}
  
