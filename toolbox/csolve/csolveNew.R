csolve <- function(FUN,x,...,gradfun=NULL,crit=1e-7,itmax=20,verbose=TRUE,alpha=1e-3,delta=1e-6){
### FUN:      Either an expression vector consisting of the equations to be solved, or a function
###           written so that if presented with a matrix x, it produces a return value of
###           same dimension as x.  In this case the number of rows in x and FUN(x) are always the
###           same.  The number of columns is the number of different input arguments
###           at which FUN is to be evaluated.  (If FUN is given as expressions, analytic derivatives are used.)
###           **************************************
### x:        initial value for FUN's argument.  A matrix with one column if FUN is a function,
###           a list setting variable names=values if FUN is an expression vector.
### gradfun:  string naming the function called to evaluate the gradient matrix.  If this
###           is NULL (the default) and FUN is a function, a numerical gradient is used instead.
###           If FUN is an expression list,
### crit:     if the sum of absolute values that FUN returns is less than this,
###           the equation is solved.
### itmax:    the solver stops when this number of iterations is reached, with rc=4
### ...:      in this position the user can place any number of additional arguments, all
###           of which are passed on to FUN and gradfun (when it is non-empty) as a list of 
###           arguments following x.
### -------------------------------------------------------------------
###           Arguments below this usually can be left at defaults.
### verbose:  If set to FALSE, the amount of output is reduced.
### alpha:    Tolerance on rate of descent.  The algorithm may produce search directions nearly
###           orthogonal to the gradient, and hence nearly zero directional derivative.  Smaller
###           alpha allows closer approach to orthogonality.
### delta:    difference interval in (cheap, forward-difference) numerical derivative.  Ignored if gradfun non-NULL.
### rc:       0 means normal solution, 1 and 3 mean no solution despite extremely fine adjustments
###           in step length (very likely a numerical problem, or a discontinuity). 4 means itmax
###           termination.
###------------ analyticg --------------
  analyticg <- !is.null(gradfun)        #if the grad argument is NULL, numerical derivatives are used.
###-------------------------------------
  nv <- length(x)
  tvec <- delta*diag(nv)
  done <- FALSE
  f0 <- FUN(x,...)
  af0 <- sum(abs(f0))
  af00 <- af0
  itct <- 0
  while(!done) {
    if((itct%%2)==1  && af00-af0 < crit*max(1,af0) && itct>3) {
      randomize <- TRUE
    } else {
      if(!analyticg) {
        grad <- (FUN(matrix(x,nv,nv)+tvec,...)-matrix(f0,nv,nv))/delta
      } else {                          # use analytic gradient
        grad <- gradfun(x,...)
      }
      if(is.real(grad)&& is.finite(grad)) {
        svdg <- svd(grad)
        if(max(svdg$d)/min(svdg$d)>1e14){
          svdg$d <- pmax(svdg$d,max(svdg$d)*1e-13)
          grad <- svdg$u%*% diag(svdg$d,ncol=length(svdg$d)) %*% t(svdg$v)
        }
        dx0 <- -solve(grad,f0)
        randomize <- FALSE
      } else {
        if(verbose){cat("gradient imaginary or infinite\n")}
        randomize <- TRUE
      }
    }
    if(randomize) {
      if(verbose){cat("Random Search\n")}
      dx0 <- sqrt(sum(x*x))/matrix(rnorm(length(x)),length(x),1)
    }
    lambda <- 1
    lambdamin <- 1
    fmin <- f0
    xmin <- x
    afmin <- af0
    dxSize <- sqrt(sum(dx0^2))
    factor <- .6
    shrink <- TRUE
    subDone <- FALSE
    while(!subDone) {
      dx <- lambda*dx0
      f <- FUN(x+dx,...)
      af <- sum(abs(f))
      if(af < afmin) {
        afmin <- af
        fmin <- f
        lambdamin <- lambda
        xmin <- x+dx
      }
      if( ((lambda >0) && (af0-af < alpha*lambda*af0)) || ((lambda<0) && (af0-af < 0) )) {
        if(!shrink) {
          factor <- factor^.6
          shrink <- TRUE
        }
        if(abs(lambda*(1-factor))*dxSize > .1*delta) {
          lambda <- factor*lambda
        } else {
          if( (lambda > 0) && (factor==.6) ) { #i.e., we've only been shrinking
            lambda <- -.3
          } else {
            subDone <- TRUE
            if(lambda > 0) {
              if(factor==.6) {
                rc <- 2
              } else {
                rc <- 1
              }
            } else {
              rc <- 3
            }
          }
        }
      } else {
        if( (lambda >0) && (af-af0 > (1-alpha)*lambda*af0) ) {
          if(shrink) {
            factor <- factor^.6
            shrink <- FALSE
          }
          lambda <- lambda/factor
        } else {                        # good value found
          subDone <- TRUE
          rc <- 0
        }
      }
    }                                   # while ~subDone
    itct <- itct+1
    if(verbose){
      cat(paste("itct ",itct,", af ",afmin,", lambda ", lambdamin,", rc ",rc),"\n")
      cat("x: \n",formatC(xmin,width=10,digits=6),"\n")
      cat("fmin:\n",formatC(fmin,width=10,digits=6),"\n")
    }
    x <- xmin
    f0 <- fmin
    af00 <- af0
    af0 <- afmin
    if(itct >= itmax) {
      done <- TRUE
      rc <- 4
    } else {
      if(af0<crit) {
        done <- TRUE
        rc <- 0
      }
    }
  }                                     #while not done
  return(list(x=xmin,f=fmin,itcount=itct,retcode=rc))
}
