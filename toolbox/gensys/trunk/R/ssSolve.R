ssSolve <- function(ex, x0, param, crit=1e-7, itmax=20, verbose=TRUE, alpha=1e-3, delta=1e-6, long=FALSE) {
  ex2 <- sssys(ex)
  fss <- derivVec(ex2, attr(ex2, "vlist"))
  csout <- csolve(fss, x0, param=param, gradfun=fss, crit=crit, itmax=itmax, verbose=verbose, alpha=alpha, delta=delta, long=long)
  return(list(csout=csout, fss=fss))
}
