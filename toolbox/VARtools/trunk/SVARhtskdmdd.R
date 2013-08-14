SVARhtskdmdd <- function(ydata,lags,xdata=NULL, const=TRUE, A0, lmd, Tsigbrk, breaks=NULL,
                         urprior=list(lambda=5,mu=1), mnprior=list(tight=3,decay=.5),
                         vprior=list(sig=NULL,w=1), train=0,flat=FALSE,nonorm=FALSE,ic=NULL)
### This gives the posterior integrated over A+ (the right-hand side coefficients), conditional
### on A0 and lmd.
###---------------------------------------------
### ydata:        endogenous variable data matrix, including initial condition dates.
### xdata:        exogenous variable data matrix, including initial condition dates.  
### const:        Constant term is added automatically if const=TRUE.
### A0:           Contemporaneous coefficient matrix --- constant.
### lmd:          Column Vectors of log variances of structural shocks.
### Tsigbrk:      Dates at which lmd vectors change.  Last date with old lmd (not first with new).
### breaks:       breaks in the data.  The first lags data points after a break are used
###               as new initial conditions, not data points for the fit.
### lambda:       weight on the co-persistence prior dummy observation.  (5 is reasonable)
###               lambda>0 => x variables included; lambda<0 => x variables excluded;
### mnprior       see vprior() comments
### urprior:      
### train:        If non-zero, this is the point in the sample at which the
###               "training sample" ends.  Prior x likelihood to this point is weighted to
###               integrate to 1, and therefore is treated as if it were itself the prior.
###               To do a pure training sample prior, set lambda=mu=0, mnprior=NULL, vprior$w=0,
###               train>lags.  
### flat:         Even with lambda=mu=vprior$w=0, mnprior=NULL, det(Sigma)^(-(nv+1)/2) is used
###               as a "prior", unless flat=TRUE. flat=TRUE is likely not to work unless train is reasonably large.
### nonorm:       If true, use dummy observations but do not normalize posterior to make them a
###               proper prior.  Useful to duplicate results obtained by others, to use
###               dummy observations that do not imply a proper prior, or to save computing time in case only the
###               posterior on this model's parameters, not the weight on the model, is needed.  
### ic:           Initial conditions matrix for use in forming the sums of coefficients dummy observations.
###               If ic=NULL, the means of the first lags observations in ydata are used.  If !is.null(ic),
###               ic should be a single "observation" on the y's and x's that will be used as the persistent
###               values entering the sums of coefficients dummies.
###
###               Note that to enter a prior directly as dummy observations, one can treat the
###               Dummy observations as a training sample.
###
{
    if (is.null(dim(ydata)))  ydata <- matrix(ydata, ncol=1)
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=list(lambda=lambda, mu=mu)) # vp$: ydum,xdum,pbreaks
    ## -------- set lmd for prior dummies --------------
    if (!is.null(dim(lmd))) {
        lmdbar <- apply(lmd, 2, mean)
    } else {
        lmdbar <- lmd
    }
    Tsigbrk <- c(invTime(Tsigbrk, ydata), T)            #dummy obs at end
    lmd <- cbind(lmd, lmdbar)
    ##-------------------------------------------
    ## var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum), breaks=matrix(c(breaks, T, T + vp$pbreaks), ncol=1),
    ## const=FALSE, lambda=lambda, mu=mu, ic=ic) # const is FALSE in this call because ones alread put into xdata
    var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum),
        breaks=matrix(c(breaks, T, T + vp$pbreaks), ncol=1), const=FALSE, lambda=urprior$lambda,
        mu=urprior$mu, ic=ic, sigpar=list(A0,lmd,Tsigbrk))
    ##  const is FALSE in this call because ones alread put into xdata
    Tu <- dim(var$u)[1]
    if ( var$snglty > 0 ) error( var$snglty, " redundant columns in rhs matrix")
    w <- -.5 * apply(var$u^2, 2, sum)  
    if(train!=0) {
        if(train <= lags)
            {
                cat("end of training sample <= # of lags\n")  #
                    return
            }
        Tp <- train
        tbreaks <- c(breaks[breaks<train],Tp)
    } else {
        Tp <- lags
        ## because need initial conditions to form lambda/mu prior dummy obs
        tbreaks <- Tp
    }
    ytrain <- ydata[1:Tp,,drop=FALSE]
    xtrain <- xdata[1:Tp,,drop=FALSE]
    priorTsigbrk <- intersect(Tsigbrk, c(1:Tp, T))
    priornsig <- length(Tsigbrk) + 1
    priorlmd <- cbind(lmd[ , 1:(priornsig - 1)], lmd[ , dim(lmd)[2]])
    if (!nonorm) { 
        varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                       breaks=c(tbreaks, Tp+vp$pbreaks), 
                       lambda=NULL, mu=NULL, const=FALSE, ic=ic,
                       sigpar=list(A0=A0,lmd=priorlmd, Tsigbrk=priorTsigbrk))
        ## const is FALSE here because xdata already has a column of ones.
        if (varp$snglty > 0) {
            warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
        } else {
            Tup <- dim(varp$u)[1]
            wp <- matrictint(crossprod(varp$u),varp$xxi,Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi)
            w=w-wp
        }
    } else {
        varp <- NULL
    }
    return(list(w=w,var=var,varp=varp,prior=list(lambda=lambda,mu=mu,vprior=vprior,mnprior=mnprior)))
}