DiscFP <- function(param, gy, y, U, alph, nit=10, crit=1e-7, theta) {
  screp <- crit + 10
  itct = 0
  while (screp > crit && itct < nit) {
    itct <- itct + 1
    Dout <- DiscPObjXnoD(param, gy, y, U=U, alph=alph, theta=theta)
    ny <- length(y)
    np <- (length(param) + 1)/2
    pold <- param[1:(np-1)]
    pold <- c(pold, 1 - sum(pold))
    pnew <- Dout$pnew
    pinew <- Dout$ygivenx %*% y * theta / (theta - 1)
    if ( length(pnew) != length(pold)) print(itct)
    screp <- sum(abs(c(pnew - pold, pinew - param[np - 1 + (1:np)])))
    param <- c(pnew[-np], pinew)
  }
  return(c(Dout=Dout, pi=pinew, itct=itct, screp=screp))
}