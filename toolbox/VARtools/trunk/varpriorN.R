varpriorN <-  function(nv=1,nx=0,lags=1,mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1),
                       urprior=list(lambda=NULL, mu=NULL), xsig=NULL, ybar, xbar=1, nstat=rep(TRUE,nv)) {
  ## varprior produces dommy observations, interpreted as scaled relative to residual variances.
  ## This function uses the var prior output to produce a mean and covariance matrix for a normal prior.
  prior <- varprior(nv, nx, lags, mnprior,vprior, urprior, xsig, ybar, xbar, nstat)
  ##------------------------------------------------------------------------
   ## contruct prior mean and variance from varprior output
  vpout <- rfvar3(prior$ydum, lags=lags, xdata=prior$xdum, const=FALSE,
                  breaks=prior$pbreaks, lambda=NULL, mu=NULL)
  shat <- with(vpout, c(rbind(matrix(aperm(By, c(2,3,1)), prod(dim(By)[2:3]), dim(By)[1]), t(Bx))))
  ## lags,vbls for y, then x, then eq
  ## Indexing in shat: ((vbl x lag), const) x eqn
  sighat <- with(vpout, kronecker(crossprod(u), xxi))
  sighat <- wtvec * t(wtvec * sighat)
  return(list(shat=shat, sighat=sighat))
}
  
  
