DiscPObj <- function(param, gx, y, U, alph) {
  ## "p" and "x" index the variables wrt which we are taking derivs
  ## "weight" indexes the nx rows of the joint density
  ## "y" indexes the ny columns of the joint density
  require(tensorA)
  nx <- (length(param)+1)/2
  p <- param[1:(nx - 1)]
  p <- to.tensor(c(p, 1 - sum(p)), c(1, "weight"))
  x <- to.tensor(param[nx:(2 * nx - 1)]. c(1, "weight"))
  ny <- length(gx)
  if ( any(p < 0) ) return(1e20)
  mm <- function(z) { matrix(z, nx, ny, byrow=TRUE) }
  sy <- function(z) {apply( z, MAR=1, FUN=sum)}
  Umat <- to.tensor(outer(x, y, FUN=U),c(weight=nx, y=ny))
  ## peaU <- p %*% exp(alph * Umat)
  eaU <- exp(alph * Umat),
  DPpeaU <- eaU[[weight=~p]]            #renames one dim
  peaU <- p %*% eaU                     #warning:  this, like some others below, is really only the
                                        #  diagonal of DPpeaU, truly indexed
                                        #  by both p and weight
  DXU <- to.tensor(attr(Umat,"gradient"), c(x=nx, y=ny))
  ##assumes U fcn written to return gradient attribute,
  ##e.g. via deriv().
  DXpeaU <- mm(p %*% (alph * eaU)) * alph * DXU
  h <-  gx / peaU
  DPh <- -DPpeaU * mm(c(gx / peaU^2))
  DXh <- -DXpeaU * mm(c(gx / peaU^2))
  eaUh <- eaU * mm(c(h))
  f <- c(p) * eaUh
  DPf <- eaUh +  p * eaU * DPh
  DXf <- c(p) * eaU * DXh + c(p) * (alph * DXU * eaUh)
  pnew <- sy(f)
  DPpnew <- sy(DPf)
  DXpnew <- sy(DXf)
  ipplus <- p > 0
  fplus <- f[ipplus, ]
  pplus <- pnew[ipplus]
  roweight <- to.tensor(sy(eaUh), c(weight=nx,1))
  ygivenx <- c(1/roweight) * eaUh
  DProweight <- to.tensor(eaU %*% t(DPh), c(weight=nx, p=nx)) #nx x nx matrix
  DXroweight <- to.tensor(sy( eaU * DXh) + (eaU * alph * DXU) %*% h, c(weight=nx, x=nx)) 
  ## result of line below is weight=nx by p=nx by y=ny tensor.  Weight indexes rows of ygivenx, p indexes what deriv's are wrt.
  DPygivenx <- -mul.tensor(DProweight,NULL , 1/roweight^2, NULL, by="weight")
               + mul.tensor(to.tensor( (1/roweight) * eaU , c(weight=nx, y=ny)), NULL, to.tensor(DPh, c(p=nx, y=ny), by="y"))
                     
  DXygivenx <- -mul.tensor(DXroweight, NULL, 1/roweight^2, NULL, by="weight")
               + mul.tensor(to.tensor( (1/rowweight) * eaU, c(weight=nx, y=ny)), NULL, to.tensor(DXh, c(p=nx, y=ny), by="y"))
    obj <- sum(f  * Umat ) - (1/alph) * pplus %*% sy(log(ygivenx[ipplus]) * ygivenx[ipplus])
  DPobj <- vector("numeric", nx)
  DXobj <- DPobj
  DPobj[ipplus] <- sy(DPf * Umat)[ipplus] - (1/alph) * 
  DXobj[ipplus] <- apply(DXU * f, MAR=1, FUN=sum)[ipplus]
  + apply(Umat * DXf, MAR=1, FUN=sum)[ipplus] -(1/alph) * apply(log(ygivenx[ipplus]) * DXf[ipplus, ], FUN=sum, MAR=1)
  obj <- -obj                           #as input to minimizer
  attr(obj, "gradient") <- -c(cbind(diag(nx - 1), rep(-1, nx - 1)) %*% DPobj, DXobj)
                                        # browser()
                                        #attr(obj, "pnew") <- pnew             #allows checking that pnew = p (as it should at optimum)
  return(obj)
}
