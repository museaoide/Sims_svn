setUpScrep <- function(crange=c(0.,0.5), nc,  pw, ww, U=function(x,y) {log(ifelse(x>0,x,1e-100)) + log(ifelse(y > x, y - x, 1e-100))}, lambda=1,
                     DxU=function(x,y) {1/x - 1/(y-x)}, lbp=1e-7, PNLTY=1) {
  ## ww,pw: Points of support and marginal probabilities for w distribution.  dw:wrange::dc:crange
  ## U:     objective function U(x,w). 
  ## DxU:   derivative of U w.r.t. its first argument. If this is NULL, numerical gradient is used.
  ## crange: range of values for x 
  ## wrange: range of values for w 
  ## lambda: multiplier on the information constraint, utility cost of capacity
  ## lbp:    lower bound on elements of dc, pc, dw, pw and 1 - upper bound on their sums.
  ## PNLTY:  weight on penalty function keeping solution away from boundary.
  if (any(diff(ww) <= 0) ) stop("ww must be increasing")
  if (any(pw < 0) || sum(pw) >1) stop("pw must be a probability vector")
  nw <- length(ww)
  ## cderiv <- !is.null(DxU)
  screpDC <- function(dc=rep(.1,9), pc=rep(.111,8)) {
    ## computes objective function E[U] minus lambda * (information constraint) - penalty for closeness to boundaries
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * crange.
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted (=1 - sum of others).
    ##
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e20))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e20))
    }
    cc <- cumsum(dc)
    pc <- c(pc, 1 - sum(pc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    cltw <- outer(cc, ww, function(a,b){a<b} )
    Umat <- outer(cc,ww, FUN=U)
    pdf <- exp(Umat)
    pdf <- pc * pdf
    h <- rep(1,nc) %*% pdf
    h <- pw / h
    pdf <- t(c(h) * t(pdf ))
    EU <- sum(pdf * Umat)
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    zerop <- pdf[cltw & (pdf < lbp)]
    zerodc <- dc[dc < lbp]
    if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc)) 
    obf <- EU - lambda *k - PNLTY * (sum((lbp / zerop - 1)^4 ) + sum((lbp / zerodc -1)^4))
    ##---------------------------------
    return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }

  dobjdpdc <- function(dc=rep(.1,9), pc=rep(.111,8)) {
    ## returns derivative and function value
    cc <- cumsum(dc)
    pc <- c(pc, 1 - sum(pc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    cltw <- outer(cc, ww, function(a,b){a<b} )
    Umat <- outer(cc,ww, U)
    eU <- exp(Umat) * cltw # basically redundant, since U=1e-100 where cltw==0
    peU <- pc %*% eU
    h <- pw/peU
    dhdp <- -t(eU %*% diag(c(h/peU))) #so w indexes rows, c indexes columns on this
    eUU <- eU * Umat
    pdf <- diag(c(pc)) %*% eU
    h <- rep(1,nc) %*% pdf
    h <- c(pw / h)
    pdf <- t(c(h) * t(pdf ))
    EU <- sum(pdf * Umat)
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    ## zerop <- pdf[cltw & (pdf < lbp)] # don't need to constrain pdf.  pc is only source of zeros
    zerop <- pc[pc < lbp]
    zerodc <- dc[dc < lbp]
    obf <- EU - lambda *k - PNLTY * (sum((lbp / zerop - 1)^4 ) + sum((lbp / zerodc -1)^4))  
    dobjdpv <- eUU %*% h + t(pc %*% eUU %*% dhdp)
    dobjdpv <-  dobjdpv - lambda * (diag(log(pc)) %*% eU %*% h + eUU %*% h + eU %*% (h * log(h)) + eU %*% h
                                    + t((pc * log(pc)) %*% eU %*% dhdp + pc %*% eUU %*% dhdp + pc %*% eU %*% diag(log(h)) %*% dhdp)
                                    + t(pc %*% eU %*% dhdp) -1 - log(pc))
    dobjdpv <- dobjdpv - PNLTY *  (pc < lbp) * (-4) * (lbp/p - 1)^3 * lbp / p^2   
    DcU <- outer(cc,ww,DxU)
    pdeuh <- diag(pc) %*% (DcU * eU) %*% diag(h)
    peuh <- diag(pc) %*% eU %*% diag(h)
    dobjdcv <- diag(pc) %*% DcU %*% h - lambda * ( apply(log(peuh) * pdeuh, 1, sum)  - pdeuh %*% log(apply(peuh,2,sum)) )
    dobjdcv <- dobjdcv - PNLTY * (dc < lbp) * (-4) * (lbp/dc - 1)^3 * lbp / dc^2
    return(list(dp=dobjdpv, dc=dobjdcv, obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }
  screpDCS <- function(x) {
    dc <- x[1:nc]
    pc <- x[-(1:nc)]
    return(-screpDC(dc,pc)$obf)
  }
  DscrepDCS <- function(x) {
    dc <- x[1:nc]
    pc <- x[-(1:nc)]
    return(c(-dobjdpdc$dc, dobjdpdc$dp[-nc]))
  }
  return(list(screpDC=screpDC, DscrepDC=dobjdpdc, screpDCS=screpDCS, DscrepDCS=DscrepDCS))
}
