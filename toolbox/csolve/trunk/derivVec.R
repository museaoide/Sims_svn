derivVec <- function(ex, x=attr(ex,"vlist"), param=vector("numeric",0), xchk=function(z){TRUE}) {
  ##    ex:    a vector of expressions, or an eqsys object
  ##     x:    a character vector of variable names w.r.t. which the expressions will be differentiated
  ##           If ex is an eqsys object, the default is for x to be the vlist attribute of ex.  
  ## param:    a named vector of parameter values that will be constant in repeated calls to fret (the returned function)
  ## xchk:     Range check function.  It returns FALSE if given an x vector that is outside the domain of definition of ex.
  ##-------------------
  ##  fret:    returned value; a function that when evaluated at a numerical x, 
  ##           returns the vector of expression values, but also, as attr(value, "grad"), the gradient matrix
  ##           It uses the fixed param values set in the call to derivVec.
  nq <- length(ex)
  nv <- length(x)
  param <- param #to put param into this namespace, so it does not need to be in every call.
  xchk <- xchk
  outf <- vector("expression",nq)
  for (iq in 1:nq) {
    f <- deriv(ex[iq], x)
    outf[iq] <- f
  }
  fret <- function(z) {
    fval <- vector("numeric", nq)
    gval <- matrix(0, nq, nv)
      zv <- as.vector(z)
      names(zv) <- if (is.null(names(z))) dimnames(z)[[1]] else names(z)
    if (xchk(zv) ) {
      for (iq in 1:nq) {
        fg <- eval(outf[iq], as.list(c(zv, param)))
        fval[iq] <- fg
        gval[iq, ] <- attr(fg, "grad")
      }
      fval <- c(fval)
      attr(fval, "grad") <- gval
    } else {                            #bad z input
      fval[] <- 1e20
      gval[] <- diag(-1, nq, nv)
      attr(fval, "grad") <- gval
    }
    names(fval) <- names(ex)
    dimnames(gval) <- list(eq=names(ex), vbl=names(zv))
    return(fval)
  }
  return(fret)
}
