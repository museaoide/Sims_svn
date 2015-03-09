Ucapm <- function(a,b) {
  gamma <- .0
  xy <- crossprod(a,b)
  u <- xy - .5 * xy^2 - gamma * crossprod(a)
}