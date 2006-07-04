\name{qz}
\alias{qz}
\title{ QZ or generalized Schur decomposition }
\description{
  Calculates the QZ decomposition of a pair of square matrices.
}
\usage{
qz(a = array(1, dim = 1, 1), b = array(1, dim = 1, 1), div = 1 + 1e-14)
}
\arguments{
  \item{a}{ left hand side input matrix }
  \item{b}{ right hand side input matrix }
  \item{div}{roots greater than div in absolute value are suppressed in
    the solution }
}
\details{
 The return code values are
 \describe{\item{0:}{all went well}
  \item{\code{j > 0}:}{the algorithm failed, but
    \code{gev[k,]} is accurate for \code{k > j}} 
  \item{\eqn{N+1}:}{some other sort of failure (\eqn{N} is the
    dimension of the matrices)}
  \item{-18:}{workspace overflow; see comments in the source if this
    occurs.}
}
}
\value{
  
  If \eqn{A} and \eqn{B} are two \eqn{N x N}  matrices, we can find
  \eqn{a,b,q,z} with \eqn{A = q a z'}, \eqn{B = q b z'}, \eqn{a,b} both
  upper triangular, and \eqn{q'q = z'z = I}, where in these formulas the
  \sQuote{'} indicates both transposition and complex conjugation.

  If \eqn{A,B} are the inputs to this function, the output is a list whose
  first four elements are \eqn{a,b,q,z} from these formulas.  In
  addition the list includes
  
  \item{rc}{return code; see Details for its interpretation}
\item{gev}{ \code{N x 2} matrix of generalized eigenvalue pairs.}
}

\references{Golub, Gene H. and Charles F. Van Loan, \emph{Matrix
    Computations}, 3d edition, Johns Hopkins University Press 1996,
  section 7.7}
\author{Christopher Sims, \email{sims@princeton.edu}, \url{http://www.princeton.edu/~sims }
\examples{
qz(matrix(c(1,0,0, 1,0,0, 1,1,1),3,3),matrix(c(0,1,0, 1,1,0, 1,1,0),3,3))
}
\keyword{math }% at least one, from doc/KEYWORDS
