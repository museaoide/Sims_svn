read.eqsys <- function(con=stdin()) {
  ix <- 0
  data <- readLines(con)
  neq <- length(data)/2 - 3
  name <- vector("character", neq)
  forward <- vector("logical", neq)
  eq <- vector("expression", neq)
  vlist <- vector("character")
  param <- vector("character")
  shock <- vector("character")
  data <- matrix(data, 2, neq + 3)
  for (iq in 1:neq) {
    nf <- scan(textConnection(data[1,iq]), what="char", quiet=TRUE)
    if (length(nf) != 2) stop(paste("equation", iq, "misformatted"))
    name[iq] <- nf[1]
    forward[iq] <- identical(nf[2],"forward")
    eq[iq] <- parse(text=data[2,iq])
  }
  vlist <- scan(textConnection(data[2, neq+1]), what="char", quiet=TRUE)
  param <- scan(textConnection(data[2, neq+2]), what="char", quiet=TRUE)
  shock <- scan(textConnection(data[2, neq+3]), what="char", quiet=TRUE)
  class(eq) <- c("eqsys", "expression")
  names(eq) <- name
  eq <- structure(eq, forward=forward, vlist=vlist, param=param, shock=shock)
  return(eq)
}
