#' Refit a model based on collapsed factors
#'
#' Refit a model based on collapsed factors
#' @param varia.list a vector of names of factors to be collapsed
#' @param merge.list a list of graycodes, each corresponding to factor in the varia.list
#' @param mod the model to be refitted
#' @return A list with the elements
#' \item{BIC}{The BIC value of the refitted model}
#' \item{AIC}{The AIC value of the refitted model}
#' \item{logLik}{The log-likelihood of the refitted model}
#' \item{refitted model}{Details of the refitted model}
#' @examples
#' fc.model.refit(varia.list=c("Kilometres"),
#'                merge.list=list(c(1,2,3,4,5)),
#'                mod = mod1.1)
#' fc.model.refit(varia.list=c("Kilometres, Make"),
#'                merge.list=list(c(1,2,3,4,5), c(1,2,3,4,5,6,7,8,8)),
#'                mod = mod1.1)
#' @export fc.model.refit
#' @export fc.model.refit.glm
#' @export fc.model.refit.rxGlm
#' @export fc.model.refit.lm

fc.model.refit <- function(varia.list,
                           merge.list,
                           mod){
  UseMethod("fc.model.refit", mod)
}

#' @describeIn fc.model.refit Refit the standard glm model class
fc.model.refit.glm <- function(varia.list,
                               merge.list,
                               mod){
  ndata <- mod$data
  varnum <- length(varia.list)
  for (each in c(1:varnum)){
    varia <- ndata[ , which(colnames(ndata) == varia.list[each])]
    new.varia <- factor.level.collapsing(varia, merge.list[[each]] )
    ndata[ , which(colnames(ndata) == varia.list[each])] <- new.varia
    if (nlevels(new.varia) == 1){
      ndata[ , which(colnames(ndata) == varia.list[each])] <- rep(1, length(new.varia))
    }
  }
  model.call <- mod$call
  model.call$data <- ndata
  glm.eval <- eval.parent(model.call)
  return(list(BIC=BIC(glm.eval), AIC=AIC(glm.eval), loglik=logLik(glm.eval), glm.eval))
}

#' @describeIn fc.model.refit Refit the rxGlm model class
fc.model.refit.rxGlm <-function(varia.list,
                               merge.list,
                               mod){

  ndata <- eval(getCall(mod)$data,
                envir = environment(formula(mod)))

  varnum <- length(varia.list)
  for (each in c(1:varnum)){
    varia <- ndata[ , which(colnames(ndata) == varia.list[each])]
    new.varia <- factor.level.collapsing(varia, merge.list[[each]] )
    ndata[ , which(colnames(ndata) == varia.list[each])] <- new.varia
    if (nlevels(new.varia) == 1){
      ndata[ , which(colnames(ndata) == varia.list[each])] <- rep(1, length(new.varia))
    }
  }
  model.call <- mod$call
  model.call$data <- ndata
  glm.eval <- eval.parent(model.call)
  aic.val <- glm.eval$aic
  bic.val <- rxAICtoBIC(glm.eval, ndata)[1]
  loglik.val <- as.numeric((glm.eval$aic - 2*glm.eval$df[1])/(-2))
  return(list(BIC=bic.val, AIC=aic.val, loglik=loglik.val,  glm.eval))
}

#' @describeIn fc.model.refit Refit the lm model class
fc.model.refit.lm <- function(varia.list,
                              merge.list,
                              mod){
  ndata <- mod$model
  varnum <- length(varia.list)
  for (each in c(1:varnum)){
    varia <- ndata[ , which(colnames(ndata) == varia.list[each])]
    new.varia <- factor.level.collapsing(varia, merge.list[[each]] )
    ndata[ , which(colnames(ndata) == varia.list[each])] <- new.varia
    if (nlevels(new.varia) == 1){
      ndata[ , which(colnames(ndata) == varia.list[each])] <- rep(1, length(new.varia))
    }
  }
  model.call <- mod$call
  model.call$data <- ndata
  lm.eval <- eval.parent(model.call)
  return(list(BIC=BIC(lm.eval), AIC=AIC(lm.eval), loglik=logLik(lm.eval), lm.eval))
}
