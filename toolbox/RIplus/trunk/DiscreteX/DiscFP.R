#' Rational inattention discrete fixed point solution
#'
#' Iterates on the first order conditions, taking the number of points of
#' support as given.
#'
#' @param param First argument of \code{DiscPObjXmv()}.  First is element is
#'           \code{nx}.  Next \code{nx-1} are marginal probabilities on 1st
#'           \code{n-1} \code{x} points. Remainder fill the
#'           \code{nx} x \code{mx} matrix of positions of the discrete x points
#' @param gy The \code{my}-length vector of  prior probabilities on \code{y}.
#' @param y  The \code{my} x \code{ny} matrix of values of \code{y}.
#' @param U  A function of \code{x,y}, where \code{x} is a length-\code{nx}
#'           vector and \code{y} is a length-ny vector, which returns the
#'           realized utility
#' @param xfcn A function of \code{ygivenx} and \code{y} that returns the
#'           choice of \code{x} for the specified distribution of \code{y}.
#' @param alph The inverse of the untility-cost of information (per nat).
#' @param nit Maximum number of iterations.
#' @param crit When the fixed point iteration produces a change smaller than
#'            this, it stops.
#' @param ... optional arguments that can be passed to \code{U()}.
#' @export 
#' ---------------------------------------
DiscFP <- function(param, gy, y, U, xfcn, alph, nit=10, crit=1e-7, ...) {
  ##   ## xfcn:   A function of ygivenx (cond'l prob matrix), and y returning optimal x
  screp <- crit + 10
  itct = 0
  while (screp > crit && itct < nit) {
    itct <- itct + 1
    ##Dout <- DiscPObjXnoD(param, gy, y, U=U, alph=alph, theta=theta)
    Dout <- DiscPObjXmv(param, gy, y, U=U, alph=alph,...)
    nx <- param[1]
    ny <- dim(y)[2]
    pold <- param[2:nx]
    pold <- c(pold, 1 - sum(pold))
    xold <- param[nx + 1:(nx * ny)]
    pnew <- Dout$pnew
    ## xnew <- Dout$ygivenx %*% y #-------xnew solves E[Dxu %*% y = 0] This is for UmvTrack
    xnew <- xfcn(Dout$ygivenx, y)
    if ( length(pnew) != length(pold)) print(itct)
    screp <- sum(abs(c(pnew - pold, xnew - xold)))
    param <- c(nx, pnew[-nx], xnew)
    if (itct %% 100 < 2) {
      print(itct)
      print(screp)
      print(Dout$obj, digits=10)
      print(Dout$pnew)
      print(xnew)
    }
  }
  return(list(Dout=Dout, x=xnew, itct=itct, screp=screp, call=match.call()))
}
