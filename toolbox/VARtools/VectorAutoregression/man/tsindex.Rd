\name{tsindex}
\alias{tsindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tsindex(dates = matrix(c(1970, 1, 1970, 2), 2, 2, byrow = TRUE), tspx = c(1960, 1990.25, 4), tsd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dates}{ ~~Describe \code{dates} here~~ }
  \item{tspx}{ ~~Describe \code{tspx} here~~ }
  \item{tsd}{ ~~Describe \code{tsd} here~~ }
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
function(dates=matrix(c(1970,1,1970,2),2,2,byrow=TRUE),tspx=c(1960,1990.25,4),tsd=NULL)
  {
    ## returns the integer row indices in the time series object tsd, or in a generic ts object A with tspx=tsp(A), corresponding
    ## to the date pairs in dates.
    if(!is.null(tsd))
      {
        tspx <- tsp(tsd)
      }
    j <- (dates[,1]-tspx[1])*tspx[3]+dates[,2]
    return(j)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
