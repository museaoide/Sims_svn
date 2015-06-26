xfcn1dlq <- function(ygivenx, y) {
    nx <- dim(y)[1]
    xhat <- ygivenx %*% y
    return(xhat)
}
