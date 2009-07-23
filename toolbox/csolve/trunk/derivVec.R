derivVec <- function(ex, x) {
  ##    ex:    a vector of expressions
  ##     x:    a character vector of variable names w.r.t. which the expressions will be differentiated
  ## param:    a vector of other variable names that will be arguments to the function returned
  ##-------------------
  ##  fret:    returned value; a function that when evaluated at a numerical x, param pair,
  ##           returns the vector of expression values, but also, as attr(value, "grad"), the gradient matrix
  nq <- length(ex)
  nv <- length(x)
  outf <- vector("expression",nq)
  for (iq in 1:nq) {
    f <- deriv(ex[iq], x)
    outf[iq] <- f
  }
  fret <- function(z, p) {
    fval <- vector("numeric", nq)
    gval <- matrix(0, nq, nv)
    zv <- as.vector(z)
    names(zv) <- dimnames(z)[[1]]
    for (iq in 1:nq) {
      fg <- eval(outf[iq], as.list(c(zv, p)))
      fval[iq] <- fg
      gval[iq, ] <- attr(fg, "grad")
    }
    fval <- c(fval)
    names(fval) <- names(ex)
    dimnames(gval) <- list(names(zv), names(zv))
    attr(fval, "grad") <- gval
    return(fval)
  }
  return(fret)
}
