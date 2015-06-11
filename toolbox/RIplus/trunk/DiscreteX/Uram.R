Uram <- function(p,c) {
    theta <- 1.5
    profit <- if( p > c) -theta * log(p) + log(p - c) else -Inf
    return(profit)
}
