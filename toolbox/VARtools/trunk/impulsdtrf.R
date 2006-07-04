impulsdtrf <- function(B,smat,nstep)
### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
###                 impulse responses" it should have the property that crossprod(smat)=sigma,
###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
###                 way to get such a smat is to set smat=chol(sigma).  To get the smat
###                 corresponding to a different ordering, use
###                 smat=chol(P %*% Sigma %*% t(P)) %*% P, where P is a permutation matrix.
###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
###                 positive powers of L in A), W=A_0^(-1)'.
###                 In general, though, it is not required that W be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
  {
    ##-----debug--------
    ##browser()
    ##------------------
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
    response <- array(0,dim=c(neq,nvar,nstep));
    response[,,1] <- t(smat);
    for (it in 2:nstep)
      {
        for (ilag in 1:min(lags,it-1))
          response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
      }
    dimnames(response) <- list(dimnames(B)[[2]],dimnames(smat)[[1]],as.character(0:(nstep-1)))
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }
