#' convert from AIC to BIC, from a rxGlm model
#'
#' convert from AIC to BIC, from a rxGlm model
#' @param rxmodel a rxGlm model
#' @param rxdata has to be the data associated with rxmode
#' @return BIC and loglikelihood
#' @export rxAICtoBIC

rxAICtoBIC <- function(rxmodel, rxdata){
  if (is.na(rxmodel$aic[1])) {stop("No AIC avaialbe for the rx-model")}
  bic.val <- rxmodel$aic - 2*rxmodel$df[1] + log(nrow(rxdata))*rxmodel$df[1]
  loglike.val <- (rxmodel$aic - 2*rxmodel$df[1])/(-2)
  return(c(bic.val, loglike.val))
}
