setUpScrep <- function(crange=c(0.,0.5), nc,  pw, ww, U=function(x,y) { log(pmax(x, .Machine$double.eps*2)) + log(pmax(y-x, .Machine$double.eps*2))}, kappa=1,
                     DxU=function(x,y)  {ifelse(x > 0 & y>x, 1/x - 1/(y-x), 0)}, lbp=1e-7, PNLTY=1) {
  ## Notes 5/28/2007:  derivative is not working?  needs to be made consistent in obf part with level fcn.  pdf emerges from csminwel not respecting c<w in first column.
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
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * crange.
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted (=1 - sum of others).
    ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
    ##        let us write the joint density as pc %*% exp(alpha*U) %*% mu(w) for a pdf pc and a weighting function mu.
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100))
    }
    cc <- cumsum(dc)
    pc <- c(pc, 1 - sum(pc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    ## cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- outer(cc,ww, FUN=U)
    ## deriv---------------
    dUmatdcc <- array(0,c(nc,nw,nc))
    for (i in 1:nc) dUmatdcc[i, , i] <- DxU(cc[i],ww)
    ## /deriv---------------
    eaU <- exp(alpha * Umat)
    ## deriv----------------
    deaUdcc <- array(eaU, c(nc,nw,nc))
    deaUdcc <- deaUdcc * dUmatdcc * alpha
    deaUdalpha <- eaU * Umat
    ## /deriv---------------
    peaU <- pc * eaU
    ## deriv----------------
    dpeaUdpc <- array(0, c(nc, nw, nc))
    for (i in 1:nc) dpeaUdpc[i, , i] <- eaU[i, ]
    dpeaUdcc <- pc * deaUdcc
    dpeaUdalpha <- peaU * Umat
    ## /deriv---------------
    pwhat <- rep(1,nc) %*% peaU
    ## deriv----------------
    dpwhatdcc <- tensor(matrix(1, 1, nc), dpeaUdcc, 2,1)
    dpwhatdpc <- tensor(matrix(1, 1, nc), dpeaUdpc, 2, 1)
    dpwhatdalpha <-  rep(1,nc) %*% dpeaUdalpha
    ## /deriv---------------
    h <- pw/pwhat
    ## deriv---------------
    dhdpwhat <-  -h/pwhat
    dhdcc <- dhdpwhat * dpwhatdcc
    dhdpc <- dhdpwhat * dpwhatdpc
    dhdalpha <- dhdpwhat * dpwhatdalpha
    ## /deriv--------------
    pdf <- t(h * t(peaU))
    ## deriv--------------
    dpdfdh <- array(0, c(nc, nw, nw))
    for (j in 1:nw) dpdfdh[ ,j, j] <- peaU[,j]
    dpdfdpc <- array(0,c(nc,nw,nc))
    eaUh <- t(h * t(eaU))
    for (i in 1:nc) dpdfdpc[i, , i] <- eaUh[i, ]
    dpdfdpc <- dpdfdpc + tensor(dpdfdh, dhdpc, 3, 1)
    dpdfdcc <- array(0,c(nc, nw, nc))
    
    dpdfdalpha <- t(c(dhdalpha) * t(peaU)) + t(h * t(dpeaUdalpha))
    ## /deriv---------------
    h <- rep(1,nc) %*% pdf
    h <- pw / h
    pdf <- t(c(h) * t(pdf ))
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