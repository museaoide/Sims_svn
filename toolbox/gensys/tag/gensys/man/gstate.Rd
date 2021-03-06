% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gstate2.R
\name{gstate}
\alias{gstate}
\title{gstate}
\usage{
gstate(G1, impact, pickS = NULL, pickC = NULL, allowRedundant = FALSE)
}
\arguments{
\item{G1}{the coefficient on  lagged y from the output of gensys.m.}

\item{impact}{the coefficient on exogenous shocks from the output of gensys.m}

\item{pickS}{an optional guess at a matrix of coefficients that extracts
the state vector from y.  Must be a matrix with nrow < ncol}

\item{pickC}{an optional guess at a matrix of coefficients that extracts
the "control" vector from y.  Must be a matrix with nrow < ncol}
}
\value{
A list with elements:
\describe{
   \item{\code{pickS}}{a matrix of coefficients that extracts a state vector
         from y.  Equal to input \code{pickS} if \code{okS} is \code{TRUE}}
   \item{\code{pickC}}{a matrix of coefficients that extracts a control
         vector from y.  Equal to input \code{pickC} if \code{okC} is
         \code{TRUE}.}
   \item{\code{GS}}{the matrix of coefficients on the lagged control and
         state variables defined by output \code{pickS} and \code{pickC}.}
    \item{\code{pickP}}{a matrix of coefficients that extracts a summary of
          past state vector.  Equal to ouput \code{pickS} if \code{okPast}
          is \code{TRUE}. This state has all the information needed at t to
          predict optimally the future of the optimized system from t+1
          onward.  May be smaller than \code{pickS}.}
    \item{\code{okS}}{input \code{pickS} works as state vector}
    \item{\code{okC}}{input \code{pickC} works with output \code{pickS} as
          control vector}
    \item{\code{okPast}}{output \code{pickS} also summarizes past (but may still be
          redundant as a summary of the past).}
    \item{\code{redundant}}{input \code{pickS} is bigger than necessary, even though
          \code{okS}.  Output \code{pickS} will be smaller, unless
          \code{allowReduntant}.}
 }
The solution was in the form \code{y(t)=G1 \%*\% y(t-1)+impact \%*\% z(t)}.  Now
it's in the form \code{pick \%*\% y(t)=GS \%*\% pick \%*\% y(t-1)+ PsiS \%*\% z(t)},
where \code{pick=rbind(pickC,pickS)}.  So in the new system, control is
stacked above state.  Also, new \code{C} is \code{HS \%*\%} new \code{S} in
solution.
}
\description{
organize a dynamic system into "state" and "control" variables
}

