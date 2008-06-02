windowA <- function(x, start=NULL, end=NULL, frequency=NULL, deltat=NULL, agg=c("average", "sum")) {
  xx <- window(x, start, end)
  tspxx <- tsp(xx)
  oldfreq <- tspxx[3]
  if (abs(oldfreq - frequency) < .01 || (is.null(frequency) && is.null(deltat)) ) {
    ## warning("No time aggregation performed")
    ## return(xx)
    noagg <- TRUE
  } else { 
    if (is.null(frequency)) {
      frequency <- round(1/deltat)
    }
    nagg <- round(oldfreq / frequency)
    if (nagg < 2) {
      noagg <- TRUE
    } else {
      tdim <- if (!is.null(dim(xx))) dim(xx)[1] else length(xx)
      vdim <- if (!is.null(dim(xx))) dim(xx)[2] else 1
      ## make sure lengths are multiples of nagg, so variables don't spill over
      leftover <- tdim %% nagg
      if (leftover != 0) {
        xx <- window(xx, end=tsp(xx)[2] + (1 - leftover / nagg) / nagg, extend=TRUE)
        tdim <- tdim + nagg - leftover
      }
      xxa <- array(as.vector(xx), c(nagg, tdim %/% nagg, vdim))
      xxa <- apply(xxa, c(2,3), sum)
      if (agg[1] == "average") xxa <- xxa/nagg
      if (leftover != 0) warning("last obs incomplete")
      xx <- ts(xxa, start=tsp(xx)[1], freq= frequency)
      noagg <- FALSE
    }
  }
  if (noagg) warning("No time aggregation")
  if (length(start(xx)) != 2 ) warning("aggregated data not aligned with time unit")

  return(xx)
}
