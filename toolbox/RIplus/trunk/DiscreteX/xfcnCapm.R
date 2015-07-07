xfcnCapm <- function(ygivenx, y, murf=.03, sigsqz=.0003) {
    ## ygivenx and y are for the risky assets only.  We fill out the riskless
    ## asset at the start here.
    require(tensor)
    nx <- dim(ygivenx)[1]
    ny <- dim(y)[1]
    y <- cbind(rep(murf, ny), y)
    my <- dim(y)[2]
    ybar <- ygivenx %*% y
    y2bar <- array(0, c(nx, my, my))
    y2 <- array(0, c(ny, my, my))
    for (iy in 1:ny)
        y2[iy, , ] <- y[iy, ] %o% y[iy, ]
    y2bar <- tensor(ygivenx, y2, 2, 1)          # nx by my by my array
    x <- matrix(0, nx, my)
    for (ix in 1:nx) {
        y2bari <- solve(y2bar[ix, , ] + sigsqz * diag(c(0, rep(1, my-1))))
        ybarix <- ybar[ix, ]
        lmd <- -(1 - sum(y2bari %*%  ybarix))/sum(y2bari)
        x[ix, ] <- y2bari %*% (ybarix - lmd)
    }
    return(x[ , -1])                      #leave out riskless weight
}
