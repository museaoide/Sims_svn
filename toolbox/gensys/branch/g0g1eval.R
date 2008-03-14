g0g1eval <- function(dexpr, x, xl, nshock, experr, param) {
  ##  dexpr:     derivative expression vector, as produced by g0g1d
  ##      x:     vector of values (not names) of x
  ##     xl:     vector of values of xl
  ## nshock:     number of exogenous shocks (all assumed centered at 0)
  ## experr:     vector of numbers or names of equations with expetational errors
  ##  param:     vector of values of parameters (all other symbols) in the expressions
  ## NB x, xl, and param must all be vectors with named components, or lists.
  nex <- length(dexpr)
  nx <- length(x)
  nxl <- length(xl)
  g0 <- matrix(0, nex, length(x))
  g1 <- matrix(0, nex, length(xl))
  Psi <- matrix(0, nex, nshock)
  Pi <- matrix(0, nex, length(experr))
  for ( i in 1:nex ) {
    dvec <- attr(eval(dexpr[i], as.list(c(x=x,xl=xl,shock=rep(0,nshock),param=param))),"gradient")
    loc <- 0
    g0[i,] <- dvec[1:nx]
    loc <- nx
    g1[i,] <- dvec[loc:(loc+nxl)]
    loc <- loc + nxl
    Psi[i,] <- dvec[loc:(loc+nshock)]
  }
  Pi[experr,] <- diag(length(experr))
  return(list(g0=g0, g1=g1, Psi=Psi, Pi=Pi))
}
