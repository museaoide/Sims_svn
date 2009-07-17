derivVec <- function(ex, x, param) {
  ##    ex:    a vector of expressions
  ##     x:    a character vector of variable names w.r.t. which the expressions will be differentiated
  ## param:    a vector of other variable names that will be arguments to the function returned
  ##-------------------
  ##  fret:    returned value; a function that when evaluated at a numerical x, param pair,
  ##           returns the vector of expression values, but also, as attr(value, "grad"), the gradient matrix
  nq <- length(ex)
  nv <- length(x)
  outf <- vector("expression",nq)
  outg <- outf
  for (iq in 1:nq) {
    f <- deriv(ex[iq], x)
    outf[iq] <- f
  }
  fret <- function(z, p) {
    names(z) <- x
    names(p) <- param
    fval <- vector("numeric", nq)
    gval <- matrix(0, nq, nv)
    for (iq in 1:nq) {
      fg <- eval(outf[iq], as.list(c(z,p)))
      fval[iq] <- fg
      gval[iq, ] <- attr(fg, "grad")
    }
    fval <- c(fval)
    attr(fval, "grad") <- gval
    return(fval)
  }
  return(fret)
}
  
