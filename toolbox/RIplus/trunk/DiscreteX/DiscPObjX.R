DiscPObjX <- function(param, gx, y, U, alph) {
  ## DiscPObj, but with useful extras returned, for use in analyzing results,
  ## not in csminwel.
  nx <- (length(param)+1)/2
  p <- param[1:(nx - 1)]
  p <- c(p, 1 - sum(p))
  x <- param[nx:(2 * nx - 1)]
  ny <- length(gx)
  if ( any(p < 0) ) return(1e20)
  Umat <- outer(x, y, FUN=U)
  peaU <- p %*% exp(alph * Umat)
  h <-  gx / peaU
  f <- c(p) * t((c(h) * t(exp(alph * Umat))))
  pnew <- apply(f, MAR=1, FUN=sum)
  fplus <- f[f>0]
  pplus <- pnew[pnew > 0]
  obj <- sum(f  * Umat ) - (1/alph) * (sum(log(fplus) * fplus) - sum(log(pplus) * pplus))
  obj <- -obj                           #as input to minimizer
  # browser()
  ygivenx <- c(h) * t(exp(alph * Umat) )
  return(list(obj=obj, pnew=pnew, ygivenx=ygivenx))
}
