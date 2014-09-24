tsregPrior <- function(vlist, lags=rep(0, length(vlist)), ldv=1, scale, 
                       bbar=c(1,rep(0,length(unlist(lags)-1))), smooth, damp, vmeans, erratio) { 
    ##     Creates dummy observations to implement a prior and provides the normalizing factor
    ##     to allow the likelihood with dummy observations to be treated as a posterior pdf.
    ## 
    ## Args:
    ##    vlist:: character vector of names of the variables
    ##    lags:   how many lags of each variable.  Can be a single integer, if all variables appear
    ##            with the same number of lags, or a vector of integers, one for each variable
    ##            [Note:  Should smooth and damp also vary with iv if lags can do so?]  
    ##    ldv:    which element of vlist is appearing as a lagged dependent variable.  Default is the first. 
    ##            If there is no lagged dependent variable, ldv=0.
    ##   scale:   a vector, of the length of vlist plus 1, of reasonable guesses as to the standard deviations of 
    ##            the variables.  The last element scales the "levels" dummy that relates the constant to the
    ##            other parameters, which should be 1.0, or smaller if guesses of variable means in vmeans are
    ##            quite uncertain.
    ##   bbar:    mean of the regression coefficient vector.  Default of 1 followed by zeros is reasonable
    ##            when ldv=1 and this is a regression with a persistent dependent variable.
    ##   smooth:  The rate at which the scale of dummy observations drops as we go from low to
    ##            high frequencies.  Should usually be in (0,1).  If close to zero, prior on the 
    ##            variable's effect is stronger at low than at high frequencies.  With smooth=1
    ##            and damp=1, the prior just shrinks all coefficients toward their prior means
    ##            symmetrically and independently.
    ##   damp:    The rate at which the dummy observations increase with increasing lag length. Should
    ##            equal or exceed one if distant lags are more likely to be large than near-in lags.}
    ##   erratio: ratio of error variance from parameter uncertainty to error variance from the
    ##            residual for the first sample observation. Keeping this constant across models of
    ##            different size avoids a prior that implies big models should fit much better.
    ##   vmeans:  The first row of the sample X matrix, including lags, but not a constant.  Or,
    ##            a vector of that length with a guess for the mean of each variable in the position
    ##            of each lagged value of that variable.  Used for setting prior on the constant.
    ##
    ## Returns:   A list with components:
    ##   X:       dummy observations implementing the prior.  Note that the last column is the
    ##            "constant", which will *not* be a column of ones.  Generating the sample
    ##            X matrix using lagts() will not by itself create a constant term in the last column.
    ##            lagts() does put the variables in the same order as this function.  But you have
    ##            to tack on the constant vector with cbind().  X is an mts object, but its
    ##            dating is arbitrary.  Use window() to give it dates that will match the
    ##            start or end of your sample. Note that lagts() allows non-sequential lists
    ##            of lags.  This program allows different lag lengths, but always 1:lag[iv].
    ##   w:       what has to be added to  the integrated likelihood, computed from the sample 
    ##            including the dummy observations, to account for the fact that the likelihood
    ##            for the dummy observations doesn't integrate to one.
    ##
    nv <- length(vlist)
    nx <- if (length(lags) > 1) {
        sum(lags)
    } else {
        nv * lags
    }
    if (length(lags) == 1) {
        lags <- rep(lags, nv)
    }
    X <- matrix(0, nx + 1, nx + 1)
    nv <- length(vlist)
    js <- 0
    for (iv in 1:nv) {
        X[js + (1:lags[iv]), js + (1:lags[iv])] <- smooth^(0:(lags[iv]-1)) * ctmat(lags[iv]) %*% diag(damp^(0:(lags[iv]-1))) * scale[iv]
        js <- js + lags[iv]
    }
    X[ ,  nx+1] <- 0
    browser()
    X[nx+1, ] <- c(vmeans, 1) * scale[nv + 1]
    erratio0 <- solve(X, c(vmeans, 1))
    erratio0 <- sum(erratio0)^2
    X <- sqrt(erratio/erratio0) * X
    w <- -.5 * determinant(crossprod(X))$modulus
    return(list(X=X, w=w, scalefac=sqrt(erratio/erratio0), call=match.call()))
}
