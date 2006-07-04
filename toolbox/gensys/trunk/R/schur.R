schur <- function(a=array(1,dim=1,1),div=1+1e-14){
  if(!is.loaded(symbol.For("zgees"))){
    dyn.load("/usr/lib/liblapack.so")}
  N<-dim(a)[1]
  a <- as.complex(a)
  SDIM<-as.integer(1)
  ev <- array(data=0+0i,dim <- c(1,N))
  z <- array(data=0+0i,dim=c(N,N))
  LWORK <- as.integer(33*N)
  WORK <- vector("complex",length=LWORK)
  RWORK <- vector("numeric",length=8*N)
  INFO <- 0
  ##-------------------------------
  ## Returns a,b,q,z s.t. q %*% t(Conj(q)) = z %*% t(Conj(z)) = I, a and b upper triangular,
  ## q %*% a %*% t(Conj(z)) = A, and q %*% b %*% t(Conj(z))=B.  The diagonal of b should be
  ## non-negative and real (a normalization).  gev is a two-column array containing the generalized
  ## eigenvalues of the pair A,B.
  ##
  ## rc=0 indicates return without problems.
  ## rc=-18 indicates not enough work space
  ## was allocated.  The curren code attempts to be overgenerous in allocating workspace
  ## based on empirical testing of the lapack code with random matrices.  If rc=-18 ever now
  ## occurs in practice, the code could be modified to first test for required space.
  ## Also, if the current version's space requirements are a problem, the code can surely
  ## be modified to use less workspace, with some possible loss in efficiency.
  ## rc=-j indicates an illegal argument type for the j'th argument to the lapack routine.  This should
  ## not occur.
  ## -------------------------------------
  ## The call to dyn.load works on a SuSE 9.0 system with the distribution's installation
  ## of lapack.  If your lapack is at a different location, modify accordingly.
  ## -------------------------------------
  ## The first two V's indicate that we want q and z computed.  The quoted
  ## "N" indicates we want no sorting.  (This would require producing a
  ## different fortran routine for each possible choice of the "stake" around
  ## which we are sorting.) The first dum is where the sort criterion routine's
  ## name would go.  N is the 
  ## dimension of all the matrices.  SDIM would be the number of eigenvalues
  ## for which the sort criterion is true, but since we're not sorting it
  ## will come back 0, always.  ALPHA and BETA are the generalized eigenvalue
  ## vectors (the diagonals of the returned a and b).  WORK and RWORK
  ## are workspace, the former complex, the latter real.  The last dum is
  ## a placeholder for logical workspace needed only if a sort is done.
  ## INFO is 0 if all goes well, -j if the j'th argument has an incorrect
  ## form, j>0 if the algorithm failed, but ALPHA(k) and BETA(k) are accurate
  ## for k>j, N+1 if there was another sort of failure.
  ## ------------------------------
  out<-.Fortran("zgees","V","N","dum",N,a,N,
                SDIM,ev,z,N,WORK,LWORK,RWORK,"dum",INFO)
  cat("workspace needed:",out[[11]][1],"\n")
  return(list(a=matrix(out[[5]],N,N),
              z=matrix(out[[9]],N,N),
              ev=out[[8]],
              rc=out[[15]]))
}
     
     
     
     
  
