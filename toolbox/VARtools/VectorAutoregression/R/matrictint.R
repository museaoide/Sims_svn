"matrictint" <-
function(S,XXi,T)
###  S:  usually sample cross product matrix of LS residuals
### XXi:  inv(X'X) matrix for rhs variables
###  T:  number of observations
###  w:  log of integrated posterior for SUR or RF VAR with det(Sigma)^(-(m+1)/2) Jeffreys-like prior
###  To get the log of the integral of the likelihood for a VAR with T observations, 
###   k rhs variables in each equation, and m equations, set T=T-m-1 and subtract .5*m*(m+1)*log(2*pi).
### We are integrating the exponential of -.5*T*m*log(2*pi)-.5*(T+m+1)*log(det(Sigma))-.5*trace(Sigma\S(beta)).
{
  k<-dim(XXi)[1]
  m<-dim(S)[1]
  cx <- try(chol(XXi))
  if(inherits(cx,"try-error")) stop("XXI not p.d.")
  cs<-try(chol(S));
  if(inherits(cs,"try-error")) stop("S not p.d.")
  ##-------debug---------
  ##browser()
  ##---------------------
  w<-(-T+k+(m-1)/2)*m*.5*log(pi)-(T-k)*sum(log(diag(cs)))+m*sum(log(diag(cx)))+ggammaln(m,(T-k)/2)
  return(w)
}
