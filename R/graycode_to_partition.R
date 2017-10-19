#' Convert graycode to listing of partition
#'
#' Convert a graycode to listing of partition (equivalence class), to illustrate which elements are grouped together.
#'
#' @param x a graycode (in any format, not necessarily in canonical format)
#' @return listing equivalence class
#' @examples
#' graycode.to.partition(c(1,2,3,1,2))
#' graycode.to.partition(c(2,1,3,2,1))
#' @seealso \code{\link[partitions]{listParts}}
#' @export graycode.to.partition

# graycode.to.partition <- function(x) {
#   x <- convert.cano.graycode(x)
#   res <- split(seq_along(x), x)
#   class(res) <- c(class(res), "equivalence")
#   return(res)
# }

graycode.to.partition <- function(x) {
  y = split(seq_along(x), x)
  z = lapply(y, paste, collapse=",")
  res <- paste(c("(",paste(z, collapse = ")(" ), ")"), collapse = "")
  return(res)
}
