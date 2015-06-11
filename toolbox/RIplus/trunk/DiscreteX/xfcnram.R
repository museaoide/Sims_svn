xfcnram <- function(cgivenp, cost) {
    theta <- 1.5
    mc <- max(cost)
    nx <- dim(cgivenp)[1]
    psup <- mc * theta / (theta - 1)
    if (any( cgivenp <= 0)) {
        cmax <- matrix(cost, nrow=nx, ncol=length(cost), byrow=TRUE)
        cmax[cgivenp <= 0] <- 0
        cmax <- apply(cmax, 1, max)
    } else {
        cmax <- rep(mc, nx)
    }
    rfcn <- function(p) {
        icx <- cost < mcp  &  cgp > 0
        theta - sum(cgp[icx] * p / (p - cost[icx]))
    }
    pstar <- vector("numeric", nx)
    for (ip in 1:nx) {
        cgp <- cgivenp[ip, ]
        mcp <- cmax[ip]
        pstar[ip] <- uniroot(rfcn, c(cmax[ip], psup), extendInt="upX")$root
    }
    pstar <- pmax(pstar, cmax)       #-Inf utility if profits negative
    return(pstar)
}
