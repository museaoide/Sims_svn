varprior <-  function(nv=1,nx=0,lags=1,mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1))
### ydum, xdum:   dummy observation data that implement the prior
### breaks:       vector of points in the dummy data after which new dummy obs start
###                   Set breaks=T+matrix(c(0,breaks),ncol=1), ydata=rbind(ydata,ydum), xdum=rbind(xdata,xdum), where 
###                   actual data matrix has T rows, in preparing input for rfvar3
### nv,nx,lags: VAR dimensions
### mnprior$tight:Overall tightness of Minnesota prior
### mnprior$decay:Standard deviations of lags shrink as lag^(-decay)
### vprior$sig:   Vector of prior modes for square roots of diagonal elements of r.f. covariance matrix
### vprior$w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
###                   vprior.sig is needed
###                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
###                   Set vprior$w=0 to achieve this.
###                   mnprior and vprior.w can each be set to NULL, thereby eliminating the corresponding
###                   dummy observations.
### Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
###                   cannot be implemented entirely with simple dummy observations.  It is also usually
###                   taken to include the sum-of-coefficients and co-persistence components
###                   that are implemented directly in rfvar3.R.  The diagonal prior on v, combined
###                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
###                   prior mean generates larger prior variances for own than for cross-effects even in 
###                   this formulation, but here there is no way to shrink toward a set of unconstrained 
###                   univariate ARs.
###-----------------------
###
{
  if (!is.null(mnprior))
    {
      xdum <- if(nx > 0) {
        array(0, dim=c(lags + 1, nx, lags, nv), dimnames=list(obsno=1:(lags + 1), xvbl=1:nx, lag=1:lags, lhsy=1:nv))
      } else {
        NULL
      }
      ydum <- array(0,dim=c(lags+1,nv,lags,nv),dimnames=list(obsno=1:(lags+1),rhsy=1:nv,lag=1:lags,lhsy=1:nv))
      for (il in 1:lags)
        {
          ##-----debug---------
          ## browser()
          ##------------------
          ydum[il+1,,il,] <- il^mnprior$decay*diag(vprior$sig,nv)
        }
      ydum[1,,1,] <- diag(vprior$sig,nv)
      ydum <- mnprior$tight * ydum
      dim(ydum) <- c(lags+1,nv,lags*nv)
      ydum <- ydum[seq(lags+1,1,by=-1),,]
      xdum <- mnprior$tight*xdum
      dim(xdum) <- c(lags+1,nx,lags*nv)
      xdum <- xdum[seq(lags+1,1,by=-1),,]
      breaks <- (lags+1)*matrix(1:(nv*lags),nv*lags,1)
      lbreak <- breaks[length(breaks)]
      breaks <- breaks[-length(breaks)] #end of sample is not a "break".  Note this makes breaks NULL if nv==lags==1.
    } else {
      ydum <- NULL;
      xdum <- NULL;
      breaks <- NULL;
      lbreak <- 0;
    }
  if (!is.null(vprior) && vprior$w>0)
    {
      ydum2 <- array(0,dim=c(lags+1,nv,nv))
      xdum2 <- array(0,dim=c(lags+1,nx,nv))
      ydum2[lags+1,,] <- diag(vprior$sig,nv)*vprior$w #The vprior$w factor was missing until 11/29/06
                                        #Original idea, not implemented, was probably that w be an integer repetition count
                                        #for variance dobs.  Now it's just a scale factor for sig.
      dim(ydum2) <- c((lags+1)*nv,nv)
      dim(ydum) <- c((lags+1)*nv,lags*nv)
      ydum <- cbind(ydum,ydum2)
      dim(xdum2) <- c((lags+1)*nx,nv)
      dim(xdum) <- c((lags +1)*nx,lags*nv)
      xdum <- cbind(xdum,xdum2)
      dim(ydum) <- c(lags+1,nv,dim(ydum)[2])
      ydum <- aperm(ydum,c(1,3,2))
      dim(ydum) <- c(dim(ydum)[1]*dim(ydum)[2],nv)
      dim(xdum) <- c(lags+1,nx,dim(xdum)[2])
      xdum <- aperm(xdum,c(1,3,2))
      dim(xdum) <- c(dim(xdum)[1]*dim(xdum)[2],nx)
      if(nv>1){
        breaks <- c(breaks, (lags+1)*(0:(nv-1))+lbreak)
      }
    } else {
      if (!is.null(ydum)) { # case with mnprior non-null, but vprior null
        ydum <- aperm(ydum, c(1, 3, 2))
        dim(ydum) <- c(prod(dim(ydum)[1:2]), dim(ydum)[3])
        xdum <- aperm(xdum, c(1,3,2))
        dim(xdum) <- c(prod(dim(xdum)[1:2]), dim(xdum)[3])
      }
    }
  return(list(ydum=ydum,xdum=xdum,pbreaks=breaks))
  ## data here in the form of T by nv y, and T x nx x.  Lagged y's not put in to a rhs regression matrix, so a "breaks" vector
  ## is needed.  
  ## rfvar3 adds persistence and sum of coeffs dummy observations to data in lhs and rhs regression matrix form.
  ## The panel VAR programs, to accommodate the connection of cs to y0's, must reorganize lhs to (T*nv) by 1 form and construct
  ## corresponding much larger (but sparse) rhs X matrix.
}
