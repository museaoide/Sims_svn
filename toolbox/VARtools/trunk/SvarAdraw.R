SvarAdraw <- function(Ahat, idmat, Sighat, SigDraws) {
  n <- dim(Ahat)[1]
  pndx <- c(t(matrix(1:(n^2),n,n)))     #permutation for vec(t())
  idndx <- (1:n^2)[idmat != 0]
  ndraw <- dim(SigDraws)[3]
  IAAI <- kronecker(A,diag(n)) + kronecker(diag(n), t(A))
  IAAI <- qr(IAAI)[ , idndx]
  dA <- qr.solve(IAAI, c(matrix(SigDraws, n^2, ndraw) - matrix(SigHat, n^2, ndraw)))
  Ad < matrix(dA, ncol=ndraw) + Ahat
  Adraw <- matrix(0, n^2, ndraw)
  Adraw[idnx, ] <- dA + c(Ahat)
  Adraw <- array(Adraw, c(n, n, ndraw))
  return(Adraw)
}
    
