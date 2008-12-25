SvarAdraw <- function(Ahat, idmat, Sighat) {
  n <- dim(Ahat)[1]
  pndx <- seq(1, n^2, by=n) %% n
  
