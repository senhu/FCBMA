#' Search the optimal model region from a close-to-optimum model
#'
#' Given a close-to-optimum model, this function continuously searches the neighbour models to find the global optimum, through a complete search
#'
#' @param merge.code the given close-to-optimum graycode, list(code1, code2)
#' @param varia.list a list of variables, c(.., .., ..)
#' @param model the full model
#' @param transition.method Transition method to change one graycode to its neighbouring graycode, choose from \code{"ChangeOne"}, \code{"GroupSplit"}. The default is \code{"ChangeOne"}.
#' @param IncludeOrigin Logical; When choose neighbouring graycode, should the input graycode be included or not.
#'
#' @return the optimal partitions
#'
#' @examples
#' \donttest{
#' data("sweden")
#' m1 <- glm(Claims ~ Kilometres+Zone+Bonus+Make, offset = log(Insured),
#'           data = sweden, family = "poisson")
#'
#' code1 <- c(1,2,2,3,3)
#' code2 <- c(1,2,3,4,5,2,6,7,8)
#' assign("addon.bic.trace", NULL, envir = environment(FCBMA))
#' assign("addon.code.trace", NULL, envir = environment(FCBMA))
#' AddOnSearch(merge.code = list(code1,code2),
#'             varia.list = c("Kilometres", "Make"),
#'             model = m1,
#'             transition.method = "ChangeOne")
#' get("addon.bic.trace", envir = environment(FCBMA))
#' get("addon.code.trace", envir = environment(FCBMA))
#'
#' assign("addon.bic.trace", NULL, envir = environment(FCBMA))
#' assign("addon.code.trace", NULL, envir = environment(FCBMA))
#' AddOnSearch(merge.code = list(c(1,2,2,2,3)),
#'             varia.list = c("Kilometres"),
#'             model = m1,
#'             transition.method = "ChangeOne")
#' get("addon.bic.trace", envir = environment(FCBMA))
#' get("addon.code.trace", envir = environment(FCBMA))
#' }
#'
#' @export AddOnSearch

AddOnSearch <- function(merge.code,
                        varia.list,
                        model,
                        transition.method,
                        IncludeOrigin = FALSE){
  UseMethod("AddOnSearch", model)
}

#' @rdname AddOnSearch
#' @export AddOnSearch.glm
#' @export
AddOnSearch.glm <- function(merge.code,
                            varia.list,
                            model,
                            transition.method,
                            IncludeOrigin = FALSE){
  current.bic <- fc.model.refit(varia.list = varia.list,
                                merge.list = merge.code,
                                mod=model)[[1]]
  num <- length(varia.list)
  all.neighbour.array <- NULL
  for (i in seq_len(num)){
    code.i <- merge.code[[i]]
    all.neighbour.i <- partition_all_neighbour(code.i,
                                               method=transition.method,
                                               IncludeOrigin = IncludeOrigin)
    all.neighbour.i <- apply(all.neighbour.i, 1, list)
    #all.neighbour.i <- list(unlist(apply(all.neighbour.i, 1, list), recursive = FALSE))
    all.neighbour.array <- append(all.neighbour.array, list(all.neighbour.i))
  }
  all.neighbour.comb <- expand.grid(all.neighbour.array)
  names(all.neighbour.comb) <- NULL
  temp.split.fun <- function(onerow){lapply(onerow, unlist)}
  all.neighbour.comb <- apply(all.neighbour.comb, 1, temp.split.fun)
  tem <- sapply(all.neighbour.comb, fc.model.refit,
                varia.list=varia.list,
                mod=model)
  neighbour.bic <- unlist(tem[1,])
  if ( current.bic > min(neighbour.bic) ){
    sel.new.neighbour <- all.neighbour.comb[which.min(neighbour.bic)][[1]]
    assign("addon.code.trace", append(get("addon.code.trace", envir = environment(FCBMA)), list(sel.new.neighbour)), envir = environment(FCBMA))
    assign("addon.bic.trace", append(get("addon.bic.trace", envir = environment(FCBMA)), min(neighbour.bic)), envir = environment(FCBMA))
    cat("number of add-on iterations:", as.integer(length(get("addon.bic.trace", envir = environment(FCBMA)))), "current min BIC:", min(neighbour.bic), "\n" )
    Recall(merge.code = sel.new.neighbour,
           varia.list=varia.list,
           model=model,
           transition.method=transition.method,
           IncludeOrigin = IncludeOrigin)
  } else {
    NewFind.BIC <- current.bic
    NewFind <- merge.code
    return(list(NewFind.BIC, NewFind))
  }
}

#' @rdname AddOnSearch
#' @export AddOnSearch.rxGlm
#' @export
AddOnSearch.rxGlm <- function(merge.code,
                              varia.list,
                              model,
                              transition.method,
                              IncludeOrigin = FALSE){
  current.bic <- fc.model.refit(varia.list = varia.list,
                                merge.list = merge.code,
                                mod=model)[[1]]
  num <- length(varia.list)
  all.neighbour.array <- NULL
  for (i in seq_len(num)){
    code.i <- merge.code[[i]]
    all.neighbour.i <- partition_all_neighbour(code.i,
                                               method=transition.method,
                                               IncludeOrigin = IncludeOrigin)
    all.neighbour.i <- list(unlist(apply(all.neighbour.i, 1, list), recursive = FALSE))
    all.neighbour.array <- append(all.neighbour.array, all.neighbour.i)
  }
  all.neighbour.comb <- expand.grid(all.neighbour.array)
  names(all.neighbour.comb) <- NULL
  temp.split.fun <- function(onerow){lapply(onerow, unlist)}
  all.neighbour.comb <- apply(all.neighbour.comb, 1, temp.split.fun)
  tem <- NULL
  for (j in seq_len(length(all.neighbour.comb))){
    tem <- c(tem, fc.model.refit(varia.list = varia.list,
                                 merge.list = all.neighbour.comb[[j]],
                                 mod=model)[[1]])
  }
  neighbour.bic <- tem
  if ( current.bic > min(neighbour.bic) ){
    sel.new.neighbour <- all.neighbour.comb[which.min(neighbour.bic)][[1]]
    assign("addon.code.trace", append(get("addon.code.trace", envir = environment(FCBMA)), list(sel.new.neighbour)), envir = environment(FCBMA))
    assign("addon.bic.trace", append(get("addon.bic.trace", envir = environment(FCBMA)), min(neighbour.bic)), envir = environment(FCBMA))
    cat("number of add-on iterations:", as.integer(length(get("addon.bic.trace", envir = environment(FCBMA)))), "current min BIC:", min(neighbour.bic), "\n" )
    Recall(merge.code = sel.new.neighbour,
           varia.list=varia.list,
           model=model,
           transition.method=transition.method,
           IncludeOrigin = IncludeOrigin)
  } else {
    NewFind.BIC <- current.bic
    NewFind <- merge.code
    return(list(NewFind.BIC, NewFind))
  }
}
