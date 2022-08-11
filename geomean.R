#' Geometric mean
#'
#' @param val Vector to take the geometric mean of
#' @export
geomean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
