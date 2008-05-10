g0g1eval <- function(dexpr, x, xl=as.vector(x), shock, experr, param) {
  ##  dexpr:     derivative expression vector, as produced by g0g1d
  ##      x:     vector of values (not names) of x
  ##     xl:     vector of values of xl 
  ## shock:      exogenous shocks (all assumed centered at 0)
  ## experr:     vector of numbers or names of equations with expectational errors
  ##  param:     vector of values of parameters (all other symbols) in the expressions
  ## NB x, xl, and param must all be vectors with named components, or lists.
  if (is.null(names(x))) { names(x) <- dimnames(x)[[1]]}
  if (is.null(names(xl))) names(xl) <- paste(names(x),"l",sep="")
  if (is.null(names(shock))) names(shock) <- dimnames(shock)[[1]]
  if (is.null(names(param))) names(param) <- dimnames(param)[[1]]
  nex <- length(dexpr)
  nx <- length(x)
  nxl <- length(xl)
  nshock <- length(shock)
  g0 <- matrix(0, nex, nx)
  g1 <- matrix(0, nex, nxl)
  Psi <- matrix(0, nex, nshock)
  Pi <- matrix(0, nex, length(experr))
  for ( i in 1:nex ) {
    dvec <- attr(eval(eval(dexpr[i]), as.list(c(xl, x, shock, param))),"gradient")
    loc <- 0
    g0[i,] <- dvec[1:nx]
    loc <- nx
    g1[i,] <- dvec[(loc+1):(loc+nxl)]
    loc <- loc + nxl
    Psi[i,] <- dvec[(loc+1):(loc+nshock)]
  }
  Pi[experr,] <- diag(length(experr))
  dimnames(g0) <- list(names(dexpr),names(x))
  dimnames(g1) <- dimnames(g0)
  if(is.numeric(experr)) experr <- names(dexpr)[experr]
  dimnames(Psi) <- list(names(dexpr),names(shock))
  dimnames(Pi) <- list(names(dexpr), experr)
  return(list(g0=g0, g1=-g1, Psi=-Psi, Pi=-Pi))
}
