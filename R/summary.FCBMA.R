#' Summarizing FCBMA method
#'
#' Summary method for class "\code{FCBMA}"
#'
#' @param object An object of class "\code{FCBMA}" resulting of a call to \code{FCBMA}.
#' @param x An object of class "\code{summary.FCBMA}", usually a result of a call to \code{summary.FCBMA}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' data("sweden")
#' m1 <- glm(Claims ~ Kilometres+Zone+Bonus+Make, offset = log(Insured),
#'           data = sweden, family = "poisson")
#' summary(m1)
#'
#' m2 <- FCBMA(m1,
#'             varia.list = c("Kilometres"),
#'             method = "complete",
#'             verbose = FALSE)
#' summary(m2)
#' }
#'
#' @export summary.FCBMA
#' @export

summary.FCBMA <- function(object, ...){
  title <- paste0("'FCBMA' with variable '",
                  paste(object$variable, collapse="', '"), "' collapsed")
  varnum <- length(object$variable)

  best.model.graycode <- object$best.state
  best.model.bic <- object$best.bic
  if (varnum==1){
    best.model.partition <- object$all.partitions[1]
  } else {
    best.model.partition <- object$all.partitions[1,]
  }
  best.model.weight <- object$all.weights[1]
  # best.fitted.model <- fc.model.refit(object$variable,
  #                                     best.model.graycode,
  #                                     mod = object$model)[[4]]
  newgraycode <- NULL
  for (i in c(1:varnum)){
    mlist <- object$table[,i]
    r1 <- NULL
    for (j in c(1:length(mlist))){
      r1 <- rbind(r1, paste(c("(",paste(mlist[[j]], collapse = ","),")"),collapse = ""))
    }
    newgraycode <- cbind(newgraycode, r1)
  }
  newtable <- object$table
  newtable[,1:varnum] <- newgraycode

  # if (best.model.summary){
  #   cat("summary of the best fitted  model:", "\n")
  #   temp.sum <- summary(x$best.fitted.model)
  #   temp.sum$call <- object$model$call
  # }

  obj <- list(title = title,
              varnum = varnum,
              variable = object$variable,
              table = newtable,
              best.model.weight = best.model.weight,
              best.model.bic = best.model.bic,
              best.model.partition = best.model.partition,
              # best.fitted.model = best.fitted.model,
              modelcall = object$model$call
  )
  class(obj) <- "summary.FCBMA"
  return(obj)
}

#' @rdname summary.FCBMA
#' @export print.summary.FCBMA
#' @export

print.summary.FCBMA <- function(x, digits = getOption("digits"), ...){
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  cat("\n")
  cat("Call: ", "\n")
  print(x$modelcall)
  cat("\n")
  for (i in c(1:x$varnum)){
    cat(paste0("Variable ", i, " '", x$variable[i],"' best collapse:"), x$best.model.partition[i] ,"\n")
  }
  cat("Best fitted model BIC:", x$best.model.bic, "\n")
  cat("Best fitted model weight:", x$best.model.weight, "\n")
  cat("\n")

  if (nrow(x$table) < 5) {
    cat("summary of the best", nrow(x$table), "collapsed models:", "\n")
    print(x$table, digits = digits)
  } else {
    cat("summary of the best 5 collapsed models:", "\n")
    print(x$table[1:5,], digits = digits)
  }
  #
  invisible(x)
}
