sylvester <- function(A, C) {
  ## solves A %*% X + X %*% t(A) + C = 0 for X.  A and C must be square.
  if(!is.loaded("dtrsyl")){
    dyn.load("/usr/lib/liblapack", now=FALSE)
  }
  n <- dim(A)[1]
  hessA <- .Fortran("dgehrd", as.integer(n), as.integer(1), as.integer(n),
                    A=as.double(A), as.integer(n), tau = double(n - 1), work = double( n * n ),
                    lwork = as.integer(n*n), info=integer(1))
  QA <- .Fortran("dorghr", as.integer(n), as.integer(1), as.integer(n), Q = as.double(hessA$A), as.integer(n),
               as.double(hessA$tau), work=double(n * n), lwork = as.integer(n * n), info=integer(1))
  schurA <- .Fortran("dhseqr", "S", "V", as.integer(n), as.integer(1), as.integer(n), T=as.double(hessA$A),
                     as.integer(n), wr=double(n), wi=double(n), q=as.double(QA$Q), as.integer(n),
                     work=double(6*n), lwork=as.integer(6*n), info=integer(1))
  q <- matrix(schurA$q, n)
  Ctran <- t(q) %*% C %*% q
  sylo <- .Fortran("dtrsyl", "N", "T", as.integer(1), as.integer(n), as.integer(n), as.double(schurA$T),
                   as.integer(n), as.double(schurA$T), as.integer(n), X=as.double(Ctran), as.integer(n),
                   scale=double(1), info=integer(1))
  X <- matrix(sylo$X, n)/sylo$scale
  return( list( X = q %*% X %*% t(q), info = sylo$info ))
}
