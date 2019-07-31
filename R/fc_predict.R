#' Prediction based on FC
#'
#' Prediction based on the fitted model given the collpase/partition that selected by running the stochastic search
#'
#' @param data new dataset, such as testing data
#' @param model model, such as training model including interactions
#' @param varia.list variables to be collapsed
#' @param merge.list graycodes
#'
#' @return prediction
#' @examples
#'
#' m1 <- glm(Claims ~ Kilometres+Zone+Bonus+Make, offset = log(Insured),
#'           data = sweden, family = "poisson")
#' a <- fc.predict(data=sweden,
#'                 model = m1,
#'                 varia.list = c("Kilometres", "Zone", "Bonus", "Make"),
#'                 merge.list = list(c(1,2,3,4,4),
#'                                   c(1,2,3,4,5,5,5),
#'                                   c(1,1,1,2,3,4,5),
#'                                   c(1,2,3,4,5,6,7,7,7)) )
#' b <- predict.glm(m1, type='response')
#' plot(cbind(a,b))
#' @export fc.predict
#' @export

fc.predict <- function(data, model, varia.list, merge.list){
  UseMethod("fc.predict", model)
}

#' @rdname fc.predict
#' @export fc.predict.glm
#' @export

fc.predict.glm <- function(data,
                           model,
                           varia.list,
                           merge.list){
  new.mod <- fc.model.refit(varia.list=varia.list,
                            merge.list=merge.list,
                            mod = model)[[4]]
  newdat <- collapseData(data, varia.list, merge.list)
  N <- nrow(newdat)
  pred <- rep(0, N)
  for (q in c(1:N)){
    tryCatch({
      newvec <- newdat[q,]
      pred[q] <- stats::predict.glm(object = new.mod, newdata = newvec, type = "response")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} )
  }
  return(pred)
}

#' @rdname fc.predict
#' @export fc.predict.rxGlm
#' @export

fc.predict.rxGlm <- function(data,
                             model,
                             varia.list,
                             merge.list){

  new.mod <- fc.model.refit(varia.list=varia.list,
                            merge.list=merge.list,
                            mod = model)[[4]]
  newdat <- collapseData(data, varia.list, merge.list)
  pred <- RevoScaleR::rxPredict(new.mod, data=newdat,
                    checkFactorLevels = FALSE, reportProgress = 0)[,1]
  return(pred)
}
