pchol <- function(sig, porder) {
  chol(sig(porder,porder))[porder,porder]
}
