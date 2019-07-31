#' Split a vector at a certain location/position
#'
#' Split a vector at a certain location/position
#'
#' @param x a vector
#' @param pos position within the vector x
#' @return a list of splited vectors
#' @examples
#' x <- c(1,2,3,4,5)
#' SplitAt(x, 3)
#'
#' @export SplitAt

SplitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% (pos+1))))
