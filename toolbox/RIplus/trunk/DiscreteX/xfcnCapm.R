xfcnCapm <- function(ygivenx, y) {
    mu <- .03                           # return on riskless asset
    gamma <- .00015  # irreducible component of variance in all assets but first
    ## ygivenx and y are for the risky assets only.  We fill out the riskless
    ## asset at the start here.
    require(tensor)
    nx <- dim(ygivenx)[1]
    ny <- dim(y)[1]
    y <- cbind(rep(mu, ny), y)
    my <- dim(y)[2]
    ybar <- ygivenx %*% y
    y2bar <- array(0, c(nx, my, my))
    y2 <- array(0, c(ny, my, my))
    for (iy in 1:ny)
        y2[iy, , ] <- y[iy, ] %o% y[iy, ]
    y2bar <- tensor(ygivenx, y2, 2, 1)          # nx by my by my array
    x <- matrix(0, nx, my)
    for (ix in 1:nx) {
        y2bari <- solve(y2bar[ix, , ] + gamma * diag(c(0, rep(1, my-1))))
        ybarix <- ybar[ix, ]
        lmd <- -(1 - sum(y2bari %*%  ybarix))/sum(y2bari)
        x[ix, ] <- y2bari %*% (ybarix - lmd)
    }
    return(x[ , -1])                      #leave out riskless weight
}
