xfcn2p <- function(ygivenx, y) {
  a <- .5  #cross-effect of other price on demand
  require(tensor)
  ybar <- ygivenx %*% y
  nx <- dim(ygivenx)[1]
  ny <- dim(y)[1]
  my <- dim(y)[2]
  M <- matrix(c(1 + a^2, -2*a, -2*a, 1 + a^2), 2)/(2 - 2*a^2)
  x <- M %*% t(ybar) + 1/(2 * (1 - a))
  return(x)
}
