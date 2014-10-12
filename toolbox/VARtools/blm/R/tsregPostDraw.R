tsregPostDraw <- function(blmout, n) {
  ncf <- length(blmout$lsout$coefficients)
  s2draw <- with(blmout, postdf/rgamma(n, postdf))
}