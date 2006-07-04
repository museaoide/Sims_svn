fdraw <- function(parvec, S, df, ndraw, T) {
  U <- solve(S)
  V <- rwwish(df, df * U / T, ndraw)
  pv <- dwish(crossprod(V), df, df * U/T)
  pard <- solve(V, rnorm(length(pard)))
  return(pard)
}
## Not sure what this is supposed to be doing.  Obviously incomplete 
## (pard undefined, parvec unused)