Ucapm <- function(theta,y) {
  gamma <- .00015                       #so z variance matches marginal y variance if sd is .01
  ## mu is yield of riskless asset
  mu <- .03
  theta <- c(1 - sum(theta), theta)
  y <- c(mu,y)
  xy <- crossprod(theta, y)
  u <- xy - .5 * xy^2 - .5 * gamma * crossprod(theta[-1])
  ## gamma is interpretable as the variance of an irreducible, zero-mean source of uncertainty
  ## in all assets except the first (risk free) asset.
  ## first version of this function applied it to all assets, and forgot the .5.
}
