print.eqsys <- function(x) {
  print(noquote(paste(names(x), ifelse(attr(x,"forward"),"*",""),":   ", as.character(as.vector(x)),sep="" )))
  invisible(x)
}
