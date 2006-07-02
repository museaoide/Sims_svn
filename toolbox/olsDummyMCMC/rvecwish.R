rvecwish <- function(v,S,cs=NULL,ndraw=1)
  ### returns ndraw draws from a factored Wishart with v d.f. and scale matrix S, or alternatively t(cs)%*%cs.
  ### The return value, when ndraw>1, is a p x p x ndraw matrix Z, where p is the number of rows in S or cs, such that for
  ### each id=1:ndraw, crossprod(Z[,,id]) is Wishart.  Note that even when ndraw=1, Z is a 3d array, so in that case crossprod(Z)
  ### is *not* Wishart.  crossprod(Z[,,1]) is.
  {
    if(is.null(cs))
      {
        if (!is.matrix(S)) 
          S <- matrix(S)
        p <- nrow(S)
        if (p != ncol(S)) {
          stop(message = "S not square in rwish().\n")
        }
        CC <- chol(S)
      }
    else
      {
        if (!is.matrix(cs)) 
          cs <- matrix(cs)
        p <- nrow(cs)
        if (p != ncol(cs)) {
          stop(message = "cs not square in rvecwish().\n")
        }
        CC <- cs
      }
    if (v < p) {
        stop(message = "v is less than the dimension of S or cs in rvecwish().\n")
    }
    if (p>1){
      Z <- array(0, dim=c(p, p,ndraw))
      chistack <- sqrt(rchisq(p*ndraw, rep(v:(v - p + 1),times=ndraw)))
      dim(chistack) <- c(p,ndraw)
      for (id in 1:ndraw){
        diag(Z[,,id,drop=FALSE]) <- chistack[,id]
      }
      dim(Z) <- c(p*p,ndraw)
      pseq <- 1:(p-1)
      px <- rep(p * pseq, pseq) + unlist(lapply(pseq, seq))
      for (id in 1:ndraw){
        pseq <- 1:(p - 1)
        Z[px,id] <- rnorm(p * (p - 1)/2)
      }
      dim(Z) <- c(p,p,ndraw)
      Z <- aperm(Z,c(2,1,3))
      dim(Z) <- c(p,p*ndraw)
      Z <- crossprod(CC,Z)
      dim(Z) <- c(p,p,ndraw)
      Z <- aperm(Z,c(2,1,3))
    }else{
      Z <- rep(0,times=ndraw)
      chistack <- sqrt(rchisq(ndraw,v))
      for (id in 1:ndraw){
        Z[id] <- chistack[id]
      }
      Z <- CC*Z
    }
    return(Z)
}
