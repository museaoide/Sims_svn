gensystvc <- function(g0, g1, c0=matrix(0,dim(g0)[1],1), psi, pi, div=-1, tvc=NULL)
  {
    ##System given as
    ##        g0*y(t)=g1*y(t-1)+c0+psi*z(t)+pi*eta(t),
    ##with z an exogenous variable process and eta being endogenously determined
    ##one-step-ahead expectational errors.  Returned system is
    ##       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) + loose*eta .
    ## If z(t) is i.i.d., the term involving fmat and fwt drops out.
    ## If the solution is unique (eu[2]==1) there is no "loose" term.  Otherwise
    ## loose characterizes the dimensions along which there is non-uniqueness.
    ## If div is omitted from argument list, a div>1 is calculated.
    ## eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu=[-2,-2] for coincident zeros.
    ## tvc:  If non-null, tvc is a matrix such that tvc %*% y(t) is a vector whose
    ## i'th element is constrained to grow slower than div[i]^t.
    ## By Christopher A. Sims 2/24/2004, from earlier matlab code of same author.
    stop("This code is under construction, not yet usable")
    eu <- c(0,0)
    realsmall <- 1e-7
    fixdiv <- (div>0)
    n <- dim(g0)[1]
    nshock <- if (is.matrix(psi)) dim(psi)[2] else if (is.null(psi)) 0 else 1
    qzl <- qz(g0,g1)
    zxz <- any((abs(diag(qzl$a))<realsmall) & (abs(diag(qzl$b))<realsmall))
    if (zxz) {
      "Coincident zeros.  Indeterminacy and/or nonexistence.\n"
      eu <- c(-2,-2)
      gev <- qzl$gev
      return(list(eu=eu,gev=gev))
    }
    zeroax <- abs(diag(qzl$a)) < realsmall
    unstabx <- abs(diag(qzl$a)) < (1-realsmall)*abs(diag(qzl$b)) # near unit roots don't count
    unstabx <- (! zeroax) & unstabx
    if (! fixdiv) {
      if (! any(unstabx)){
        div <- 1.01
      } else
      {
        div <- .5*(min(abs(diag(qzl$b)[unstabx]/diag(qzl$a)[unstabx]))+1)
      }
    }
    unstabx <- div*abs(diag(qzl$a))<= abs(diag(qzl$b))
    nunstab <- sum(unstabx)
    qzl <- qzdiv(div,qzl)
    qq <- t(Conj(qzl$q))                # to match matlab convention 
    gev <- qzl$gev
    ## note that this means that gev is not simply the diagonals of a nd b.  qzdiv
    ## changes the numbers on the diagonals (though not their ratios), but merely reorders
    ## the original gev.  
    if (nunstab==n){
      six <- NULL
      uix <- 1:n
    } else
    {
      if (nunstab==0){
        uix <- NULL
        six <- 1:n
      } else
      {
        uix <- (n-nunstab+1):n
        six <- 1:(n-nunstab)
      }
    }
    q1 <- qq[six,,drop=FALSE]
    q2 <- qq[uix,,drop=FALSE]
    z1 <- t(Conj(qzl$z[,six,drop=FALSE]))
    z2 <- t(Conj(qzl$z[,uix,drop=FALSE]))
    a2 <- qzl$a[uix,uix,drop=FALSE]
    b2 <- qzl$b[uix,uix,drop=FALSE]
    ## debug
    ## browser()
    etawt <- q2 %*% pi
    neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
    ndeta <- min(nunstab,neta)
    if(ndeta==0){
      ueta <- matrix(0,nunstab,0)
      deta <- vector("numeric",0)
      veta <- matrix(0,neta,0)
      bigev <- vector("logical",0)
    } else {
      sd <- svd(etawt)
      ueta <- sd$u; deta <- sd$d; veta <- sd$v
      bigev <- deta>realsmall
      ueta<-ueta[,bigev,drop=FALSE]
      veta<-veta[,bigev,drop=FALSE]
      deta<-deta[bigev]
    }
    eu[1] <- sum(bigev) >= nunstab
    ##----------------------------------------------------
    ## Note that existence and uniqueness are not just matters of comparing
    ## numbers of roots and numbers of endogenous errors.  These counts are
    ## reported below because usually they point to the source of the problem.
    ##------------------------------------------------------
    etawt1 <- q1 %*% pi
    ndeta1 <- min(n-nunstab,neta)
    if(ndeta1==0){
      ueta1 <- matrix(0,n-nunstab,0)
      deta1 <- vector("numeric",0)
      veta1 <- matrix(0,neta,0)
      bigev1 <- vector("logical",0)
    } else {
      sd <- svd(etawt1)
      ueta1<-sd$u
      deta1 <- sd$d
      veta1 <- sd$v
      bigev1 <- deta1 > realsmall
    }
    if (any(bigev1)) { #needed because empty dimensions are dropped after select
      ueta1 <- ueta1[,bigev1,drop=FALSE]
      veta1 <- veta1[,bigev1,drop=FALSE]
      deta1 <- deta1[bigev1]
      loose <- veta1-veta %*% t(Conj(veta)) %*% veta1
      svdl <- svd(loose)
      loose <- sum(abs(svdl$d)>realsmall*n)
      unq <- (loose==0)
    } else {
      ueta1 <- matrix(1,n-nunstab,0)
      veta1 <- matrix(1,neta,0)
      deta1 <- vector("complex",0)
      unq <- TRUE
    }
    if (unq) {
      eu[2] <- 1
    } else
    {
      cat("Indeterminacy.", loose, "loose endog errors.\n")
    }
    ## Note: if v is a vector of length n and m is an nxp matrix,
    ## v*m==diag(v)%*%m, m/v==solve(diag(v),m)==diag(v)\m (matlab notation)
    ##
    tmat <- cbind(diag(n-nunstab),
                  -t(Conj((ueta %*% (t(Conj(veta))/deta)) %*% veta1 %*% (deta1 * t(Conj(ueta1)))))  )  
    G0<- rbind( tmat %*% qzl$a, cbind(matrix(0,nunstab,n-nunstab), diag(nunstab)))
    G1<- rbind(tmat %*% qzl$b, matrix(0,nunstab,n))
    ##----------------------
    ## G0 is always non-singular because by construction there are no zeros on
    ## the diagonal of a[1:(n-nunstab),1:(n-nunstab)], which forms G0's ul corner.
    ##-----------------------
    G0I <- solve(G0)
    G1 <- G0I%*%G1
    ##----------- uix can be empty, e.g. in indeterminate systems with no unstable roots ------------
    if(is.null(uix)){
      C <- G0I %*% tmat %*% qq %*% c0
      fmat <- matrix(0,0,0)
      fwt <- matrix(0, 0, nshock)
      impact <- G0I %*% tmat %*% qq %*% psi
    }else{
      C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(qzl$a[uix,uix,drop=FALSE]-qzl$b[uix,uix,drop=FALSE],q2%*%c0) )
      impact <- G0I %*% rbind(tmat %*% qq %*% psi, matrix(0,nunstab, nshock))
      fmat <- solve(qzl$b[uix,uix,drop=FALSE],qzl$a[uix,uix,drop=FALSE])
      fwt <- -solve(qzl$b[uix,uix,drop=FALSE],q2 %*% psi)
    }
    ywt <- G0I[,uix,drop=FALSE]
    loose <- etawt1 %*% (diag(neta) - veta %*% t(Conj(veta))) 
    ##-------------------- above are output for system in terms of z'y -------
    G1<-Re(qzl$z %*% G1 %*% t(Conj(qzl$z)))
    C <- Re(qzl$z%*%C)
    impact <- Re(qzl$z%*%impact)
    ywt <- qzl$z%*%ywt
    loose <- Re(qzl$z %*% rbind(loose,matrix(0,nunstab,neta)))
    return(list(G1=G1,C=C,impact=impact,fmat=fmat,fwt=fwt,ywt=ywt,gev=gev,eu=eu,loose=loose))
  }

