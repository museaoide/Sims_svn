c.eqsys <- function(...,recursive=FALSE) {
  structure(NextMethod(...), class=c("eqsys","expression"))
}
