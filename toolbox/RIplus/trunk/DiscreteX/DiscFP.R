DiscFP <- function(param, gy, y, U, xfcn, alph, nit=10, crit=1e-7, ...) {
  ##   ## xfcn:   A function of ygivenx (cond'l prob matrix), y and U returning optimal x
  screp <- crit + 10
  itct = 0
  while (screp > crit && itct < nit) {
    itct <- itct + 1
    ##Dout <- DiscPObjXnoD(param, gy, y, U=U, alph=alph, theta=theta)
    Dout <- DiscPObjXmv(param, gy, y, U=U, alph=alph,...)
    ny <- length(y)
    nx <- param[1]
    pold <- param[2:nx]
    pold <- c(pold, 1 - sum(pold))
    pnew <- Dout$pnew
    ## xnew <- Dout$ygivenx %*% y #-------xnew solves E[Dxu %*% y = 0] This is for UmvTrack
    xnew <- xfcn(Dout$ygivenx, y, U)
    if ( length(pnew) != length(pold)) print(itct)
    screp <- sum(abs(c(pnew - pold, xnew - param[-(1:nx)])))
    param <- c(nx, pnew[-nx], xnew)
  }
  return(list(Dout=Dout, x=xnew, itct=itct, screp=screp))
}