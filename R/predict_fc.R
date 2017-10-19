#' Prediction based on FC
#'
#' Prediction based on the fitted model given the collpase/partition that selected by running the stochastic search
#' @param data new dataset, such as testing data
#' @param model model, such as training model including interactions
#' @param varia.list variables to be collapsed
#' @param merge.list graycodes
#' @return prediction
#' @examples
#' a <-predict.fc(data=dataset,
#'                model = mod1.1,
#'                varia.list = c("Kilometres", "Zone", "Bonus", "Make"),
#'                merge.list = list(c(1,2,3,4,4),
#'                                  c(1,2,3,4,5,5,5),
#'                                  c(1,1,1,2,3,4,5),
#'                                  c(1,2,3,4,5,6,7,7,7)) )
#' b <- predict.glm(new.mod, type='response')
#' @export predict.fc
#' @export predict.fc.glm
#' @export predict.fc.rxGlm

predict.fc <- function(data, model, varia.list, merge.list){
  UseMethod("predict.fc", model)
}

#' @describeIn predict.fc for standard GLM
predict.fc.glm <- function(data,
                           model,
                           varia.list,
                           merge.list){
  new.mod <- fc.model.refit(varia.list=varia.list,
                            merge.list=merge.list,
                            mod = model)[[4]]
  newdat <- collapse.data(data, varia.list, merge.list)
  N <- nrow(newdat)
  pred <- rep(0, N)
  for (q in c(1:N)){
    tryCatch({
      newvec <- newdat[q,]
      pred[q] <- predict(new.mod, newvec, type='response')
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} )
  }
  return(pred)
} #end of function

#' @describeIn predict.fc for class rxGlm in RevoScaleR package
predict.fc.rxGlm <- function(data,
                             model,
                             varia.list,
                             merge.list){

  new.mod <- fc.model.refit(varia.list=varia.list,
                            merge.list=merge.list,
                            mod = model)[[4]]
  newdat <- collapse.data(data, varia.list, merge.list)
  pred <- rxPredict(new.mod, data=newdat,
                    checkFactorLevels = FALSE, reportProgress = 0)[,1]
  return(pred)
}
