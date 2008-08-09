print.eqsys <- function(x) {
  print(noquote(paste(names(x),":   ", as.character(as.vector(x)),sep="" )))
  invisible(x)
}
