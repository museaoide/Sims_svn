setUpScrep <- function(crange=c(0.,0.5), nc,  pw, ww,
                       U=function(x,y) { log(pmax(x, .Machine$double.eps*2)) + log(pmax(y-x, .Machine$double.eps*2))}, kappa=1,
                     DxU=function(x,y)  {ifelse(x > 0 & y>x, 1/x - 1/(y-x), 0)}, lbp=1e-7, PNLTY=1) {
  ## Notes 5/28/2007:  derivative is not working?  needs to be made consistent in obf part with level fcn.
  ## pdf emerges from csminwel not respecting c<w in first column.
  ## Notes 6/1/2007:  Still need to fix derivative.  With numgrad seems possibly to be working, but numerically rough.
  ## ----------------FIX ABOVE---------------------
  ## ww,pw: Points of support and marginal probabilities for w distribution.  
  ## U:     objective function U(x,w). 
  ## DxU:   derivative of U w.r.t. its first argument. If this is NULL, numerical gradient is used.
  ## crange: range of values for x 
  ## nc:     number of distinct c values.  (length of dc)
  ## lambda: multiplier on the information constraint, utility cost of capacity
  ## lbp:    lower bound on elements of dc, pc, dw, pw and 1 - upper bound on their sums.
  ## PNLTY:  weight on penalty function keeping solution away from boundary.
  if (any(diff(ww) <= 0) ) stop("ww must be increasing")
  if (any(pw < 0) || sum(pw) >1) stop("pw must be a probability vector")
  nw <- length(ww)
  ## cderiv <- !is.null(DxU)
  screpDC <- function(dc=rep(.1,9), pc=rep(.111,8), alpha=1) {
    ## computes objective function, including the penalty defining the constraints, and its derivative vector
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * (crange[2] - crange[1]) + crange[1].
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted ( = 1 - sum of others).
    ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
    ##        let us write the joint density as pc %*% exp(alpha*U) %*% mu(w) for a pdf pc and a weighting function mu.
    require("tensorA")
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100))
    }
    cc <- to.tensor(cumsum(dc),dims=c(ic=nc))
    pc <- to.tensor(c(pc, 1 - sum(pc), dims=c(ic=nc)))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    ww <- to.tensor(ww,dims=c(iw=nw))
    pw <- to.tensor(pw,dim=c(iw=nw))
    ## cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- to.tensor(outer(cc, ww, FUN=U), dims=c(dim(cc), dim(ww)), what=c(1,2))
    ## deriv---------------
    dUmatdcc <- to.tensor(0,c(ic=nc, iw=nw, idc=nc))
    for (i in 1:nc) dUmatdcc[i, , i] <- DxU(cc[i],ww)
    ## /deriv---------------
    eaU <- exp(alpha * Umat)
    ## deriv----------------
    deaUdcc <- rep(eaU, nc, pos=3, name="idc")
    deaUdcc <- deaUdcc * dUmatdcc * alpha
    deaUdalpha <- eaU * Umat
    ## /deriv---------------
    peaU <- mul.tensor(X=pc, Y=eaU, by="ic") 
    ## deriv----------------
    dpeaUdpc <- to.tensor(0, c(ic=nc, iw=nw, idc=nc))
    for (i in 1:nc) dpeaUdpc[i, , i] <- eaU[i, ]
    dpeaUdcc <- mul.tensor(X=pc, Y=deaUdcc, by="ic")
    dpeaUdalpha <- peaU * Umat
    ## /deriv---------------
    pwhat <- margin.tensor(peaU, "ic")
    ## deriv----------------
    dpwhatdcc <- margin.tensor(dpeaUdcc, "ic")
    dpwhatdpc <- margin.tensor(dpeaUdpc, "ic")
    dpwhatdalpha <-  margin.tensor(dpeaUdalpha,"ic")
    ## /deriv---------------
    h <- pw/pwhat
    ## deriv---------------
    dhdpwhat <-  -h/pwhat
    dhdcc <- mul.tensor(dhdpwhat, Y=dpwhatdcc, by="iw")
    dhdpc <- mul.tensor(dhdpwhat * dpwhatdpc
    dhdalpha <- dhdpwhat * dpwhatdalpha
    ## /deriv--------------
    pdf <- mul.tensor(h, NULL, peaU, NULL, by="iw")
    ## deriv--------------
    dpdfdh <- to.tensor(0, c(ic=nc, iw=nw, idw=nw))
    for (j in 1:nw) dpdfdh[ ,j, j] <- peaU[,j]
    dpdfdpc <- to.tensor(0,c(ic=nc,iw=nw,idc=nc))
    eaUh <- mul.tensor(h, NULL, eaU, NULL, by="iw")
    for (i in 1:nc) dpdfdpc[i, , i] <- eaUh[i, ]
    dpdfdpc <- dpdfdpc + mul.tensor(dpdfdh, "idw", dhdpc, "iw")
    dpdfdcc <- to.tensor(0,c(ic=nc, iw=nw, idc=nc))
    for (i in 1:nc) dpdfdcc[i, , i] <- 
    dpdfdalpha <- mul.tensor(dhdalpha, NULL, peaU, NULL, by="iw") + mul.tensor(h, NULL, dpeaUdalpha, NULL, by="iw")
    # #t(c(dhdalpha) * t(peaU)) + t(h * t(dpeaUdalpha)) 
    ## /deriv---------------
    ##!! still need to add deriv sections for redefined h and renormed pw, plus for penalty terms.
    h <- margin.tensor(pdf, "ic")     #rep(1,nc) %*% pdf 
    h <- pw / h
    pdf <- mul.tensor(h, NULL, pdf, NULL, by="iw")   #t(c(h) * t(pdf ))
    EU <- sum(pdf * Umat)
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    zerop <- pdf[ cltw & (pdf < lbp)]
    ## zerodc <- dc[dc < lbp]
    ## if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
    zerodc <- dc[dc<0]
    if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
    ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
    ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
    obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum(zerodc^2)  + 10 * max(k-kappa,0)^2)
    ##---------------------------------
    return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }
}

##!! below is just above with deriv sections removed.
screpDO <-   function(dc=rep(.1,9), pc=rep(.111,8), alpha=1) {
  ## computes objective function only, including the penalty defining the constraints.
  ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
  ##        sum to one.  Support point i is at (sum_1^i dc_i) * (crange[2] - crange[1]) + crange[1].
  ## pc:    Marginal probabilities on the points in c's support, with last one omitted ( = 1 - sum of others).
  ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
  ##        let us write the joint density as pc %*% exp(alpha*U) %*% h(w) for a pdf pc and a weighting function h.
  require("tensorA")
  if ((sum(dc) > 1) || any( dc < 0 ) ) {
    message("bad dc vector")
    return(list(obf=-1e100))
  }
  if ((sum(pc) > 1) || any (pc < 0 ) ) {
    message("bad pc vector")
    return(list(obf=-1e100))
  }
  cc <- to.tensor(cumsum(dc),dims=c(ic=nc))
  pc <- to.tensor(c(pc, 1 - sum(pc), dims=c(ic=nc)))
  cc <- cc * (crange[2] - crange[1]) + crange[1]
  ww <- to.tensor(ww,dims=c(jw=nw))
  pw <- to.tensor(pw,dim=c(jw=nw))
  ## cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
  Umat <- to.tensor(outer(cc, ww, FUN=U), dims=c(dim(cc), dim(ww)), what=c(1,2))
  ## deriv---------------
  ##     dUmatdcc <- to.tensor(0,c(ic=nc, iw=nw, idc=nc))
  ##     for (i in 1:nc) dUmatdcc[i, , i] <- DxU(cc[i],ww)
  ## /deriv---------------
  eaU <- exp(alpha * Umat)
  ## deriv----------------
  ##     deaUdcc <- to.tensor(array(eaU, c(ic=nc, iw=nw, idc=nc)), dims=c(ic=nc, iw=nw, idc=nc), what=1:3)
  ##     deaUdcc <- deaUdcc * dUmatdcc * alpha
  ##     deaUdalpha <- eaU * Umat
  ## /deriv---------------
  peaU <- mul.tensor(pc, eaU, by="ic")
  ## deriv----------------
  ##     dpeaUdpc <- to.tensor(0, c(ic=nc, iw=nw, idc=nc))
  ##     for (i in 1:nc) dpeaUdpc[i, , i] <- eaU[i, ]
  ##     dpeaUdcc <- pc * deaUdcc
  ##     dpeaUdalpha <- peaU * Umat
  ## /deriv---------------
  pwhat <- margin.tensor(peaU, "ic")
  ## deriv----------------
  ##     dpwhatdcc <- margin.tensor(dpeaUdcc, "ic")
  ##     dpwhatdpc <- margin.tensor(dpeaUdpc, "ic")
  ##     dpwhatdalpha <-  margin.tensor(dpeaUdalpha,"ic")
  ## /deriv---------------
  h <- pw/pwhat
  ## deriv---------------
  ##     dhdpwhat <-  -h/pwhat
  ##     dhdcc <- dhdpwhat * dpwhatdcc
  ##     dhdpc <- dhdpwhat * dpwhatdpc
  ##     dhdalpha <- dhdpwhat * dpwhatdalpha
  ## /deriv--------------
  pdf <- mul.tensor(h, NULL, peaU, NULL, by="iw")
  ## deriv--------------
  ##     dpdfdh <- to.tensor(0, c(ic=nc, iw=nw, idw=nw))
  ##     for (j in 1:nw) dpdfdh[ ,j, j] <- peaU[,j]
  ##     dpdfdpc <- to.tensor(0,c(ic=nc,iw=nw,idc=nc))
  ##     eaUh <- mul.tensor(h, NULL, eaU, NULL, by="iw")
  ##     for (i in 1:nc) dpdfdpc[i, , i] <- eaUh[i, ]
  ##     dpdfdpc <- dpdfdpc + mul.tensor(dpdfdh, "idw", dhdpc, "iw")
  ##     dpdfdcc <- to.tensor(0,c(ic=nc, iw=nw, idc=nc))
  ##     dpdfdalpha <- mul.tensor(dhdalpha, NULL, peaU, NULL, by="iw") + mul.tensor(h, NULL, dpeaUdalpha, NULL, by="iw")
                                        # #t(c(dhdalpha) * t(peaU)) + t(h * t(dpeaUdalpha)) 
  ## /deriv---------------
  h <- margin.tensor(pdf, ic) #rep(1,nc) %*% pdf #Is this redundant?  Looks like should have h == 1, so pdf unchanged.
  h <- pw / h
  pdf <- mul.tensor(h, NULL, pdf, NULL, by="iw") #t(c(h) * t(pdf ))
  EU <- sum(pdf * Umat)
  k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
  ## obf <- EU - lambda * k
  ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
  ## changes in large negative values of log p's, and may speed convergence.
  ##---------------------------------
  zerop <- pdf[ cltw & (pdf < lbp)]
  ## zerodc <- dc[dc < lbp]
  ## if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
  zerodc <- dc[dc<0]
  if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
  ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
  ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
  obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum(zerodc^2)  + 10 * max(k-kappa,0)^2)
  ##---------------------------------
  return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
}
