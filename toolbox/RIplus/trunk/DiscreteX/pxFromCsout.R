pxFromCsout <- function(cx, nx) {
    p <- cx[1:(nx-1)]
    p <- c(p, 1 - sum(p))
    x <-matrix(cx[-(1:(nx-1))], nrow=nx)
    sx <- apply(x, 1, sum)
    x <- cbind(x, 1 - sx)
    return(list(p=p, x=x))
}
