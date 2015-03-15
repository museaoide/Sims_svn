xfcnCapm <- function(ygivenx, y) {
  gamma <- .001  # penalty on leverage 
  ## require(tensor)
  ybar <- ygivenx %*% y
  nx <- dim(ygivenx)[1]
  ny <- dim(y)[1]
  my <- dim(y)[2]
  y2bar <- array(0, c(nx, my, my))
  y2 <- array(0, c(ny, my, my))
  for (iy in 1:ny) y2[iy, , ] <- y[iy, ] %o% y[iy, ]
  y2bar <- tensor(ygivenx, y2, 2, 1)          # nx by my by my array
  x <- matrix(0, nx, my)
  for (ix in 1:nx) {
    y2bari <- solve(y2bar[ix, , ] + gamma * diag(my))
    ybarix <- ybar[ix,]
    lmd <- -(1 - sum(y2bari %*%  ybarix))/sum(y2bari)
    x[ix, ] <- y2bari %*% (ybar[ix,] - lmd)
  }
  return(x)
}