\name{matrictint}
\alias{matrictint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matrictint(S, XXi, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{ ~~Describe \code{S} here~~ }
  \item{XXi}{ ~~Describe \code{XXi} here~~ }
  \item{T}{ ~~Describe \code{T} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
