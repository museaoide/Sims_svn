screpDCS <- function(...){
  ## a wrapper for screpD(), to let it work as the function argument to csminwel.
  value=-screpD(...)$obf; return(ifelse(is.finite(value),value,1e20))
}
