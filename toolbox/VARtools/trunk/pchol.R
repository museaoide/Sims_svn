pchol <- function(sig, porder) {
  invporder <- match(1:length(porder), porder)
  chol(sig[porder, porder])[invporder, invporder]
}
