xfcnpolar <- function(ygivenx, y) {
    nx <- dim(ygivenx)[1]
    yrect <- cbind(y[ , 1] * cos(y[ , 2]), y[ , 1] * sin(y[ , 2]))
    yhat <- ygivenx %*% yrect
    r <- sqrt(yhat[ , 1]^2 + yhat[ , 2]^2)
    theta <- atan2(yhat[ , 2], yhat[ , 1])
    theta <- ifelse(theta < 0, 2 * pi + theta, theta)
    yopt <- cbind(r, theta)
    return(yopt)
}
    
