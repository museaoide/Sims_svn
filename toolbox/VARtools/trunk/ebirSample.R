ebirSample <- function(pdout, horiz, shockf=NULL, ... ) {
  ## pdout is the output list from postdraw(), which has
  ## elements By (equations by variables by lags by draws),
  ## Bx (equations by nx by draws), and smat (equations by
  ## equations by draws).
  ## We suppose draws are from a reduced form VAR and
  ## shockf is a one-one transformation of smat that
  ## for each i returns a square root of
  ## crossprod(t(smat[ , , i])), i.e. crossprod(t(shockf(s, ...))) ==
  ## crossprod(t(s)).  "..." is extra arguments for shockf.
  ##-------------
  ndraw <- dim(pdout$By)[4]
  if (is.null(shockf)) {
    sfac <- pdout$smat
  } else {
    sfac <- array(0, dim(pdout$smat))
    for (id in 1:ndraw) 
    sfac[ , , id] <- shockf(pdout$smat[,,id], ...)
  }
  By <- pdout$By
  nv <- dim(By)[1]
  lags <- dim(By)[3]
  ndraw <- dim(By)[4]
  resp <- array(0, c(nv, nv, horiz, ndraw))
  for (id in 1:ndraw)
    resp[ , , , id] <-   impulsdtrf(vout=list(By=mxpd$By[, , , id]), smat=mxpd$smat[, , id], nstep=horiz)
  return(resp)
}
