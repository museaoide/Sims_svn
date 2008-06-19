print.eqsys <- function(x) {
  xout <- noquote(as.character(x))
  names(xout) <- names(x)
  print(xout)
  invisible(x)
}
