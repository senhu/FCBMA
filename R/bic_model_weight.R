#' Approximate model posterior probability using BIC.
#'
#' This function allows you to calculate/approximate model posterior probability given BIC values, which can then be used for model weights in BMA weighted average.
#' @param BIC.vec a vector of all BIC values used to calculate model weight
#' @return the approximated posterior model probability, i.e. model weight used in BMA.
#' @examples
#' BIC_BMAweight(c(10,11))
#' @seealso \code{\link[BMA]{bma.glm}}
#' @export bic.model.weight

bic.model.weight <- function(BIC.vec){
  newvec <- BIC.vec-(min(BIC.vec))
  expvec <- exp(-0.5 * newvec)
  weight.vec <- expvec / sum(expvec)
  return(weight.vec)
}
