Ucapm <- function(theta,y, murf=.03, sigsqz=.0003) {
  ## mu is yield of riskless asset
  theta <- c(1 - sum(theta), theta)
  y <- c(murf,y)
  xy <- crossprod(theta, y)
  u <- xy - .5 * xy^2 - .5 * sigsqz * crossprod(theta[-1])
  ## sigsqz is interpretable as the variance of an irreducible, zero-mean source of uncertainty
  ## in all assets except the first (risk free) asset.
  ## first version of this function applied it to all assets, and forgot the .5.
}
