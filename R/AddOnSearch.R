#' Search the optimal model region from a close-to-optimum model
#'
#' Given a close-to-optimum model, this function continuously searches the neighbour models to find the global optimum, through a complete search
#' @param merge.code the given close-to-optimum graycode, list(code1, code2)
#' @param varia.list a list of variables, c(.., .., ..)
#' @param model the full model
#' @return the optimal partitions
#' @note check "addon.code.trace" and "addon.bic.trace" in global environment
#' @examples
#' code1 = c(1,2,2,3,3)
#' code2 = c(1,2,3,4,5,2,6,7,8)
#' code1 = c(1,1,1,2,3)
#' code2 = c(1,2,3,3,5,2,6,7,7)
#'
#' addon.code.trace <<- NULL
#' addon.bic.trace <<- NULL
#' AddOnSearch.function(merge.code = list(code1,code2),
#'             varia.list = c("Kilometres", "Make"),
#'             model = mod2.4)
#' addon.code.trace
#' addon.bic.trace
#'
#' addon.code.trace <<- NULL
#' addon.bic.trace <<- NULL
#' AddOnSearch.function(merge.code = list(c(1,2,2,2,3)),
#'                      varia.list = c("Kilometres"),
#'                      model = mod2.4)
#' addon.code.trace
#' addon.bic.trace
#'
#' resss <- rxAddOnSearch.function(final.merge.code = list(as.vector(Extract.SA.Res( countyfreq.res[[2]], 0.6)$State[[2]])),
#' varia.list = c("COUNTYNAME"),
#' mod=fin.freq.mod.rx)
#' @export AddOnSearch
#' @export AddOnSearch.glm
#' @export AddOnSearch.rxGlm

AddOnSearch <- function(merge.code,
                        varia.list,
                        model,
                        transition.method,
                        IncludeOrigin){
  UseMethod("AddOnSearch", model)
}

AddOnSearch.glm <- function(merge.code,
                            varia.list,
                            model,
                            transition.method,
                            IncludeOrigin){
  current.bic <- fc.model.refit(varia.list = varia.list,
                                merge.list = merge.code,
                                mod=model)[[1]]
  num <- length(varia.list)
  all.neighbour.array <- NULL
  for (i in seq_len(num)){
    code.i <- merge.code[[i]]
    all.neighbour.i <- partition.all.neighbour(code.i, method=transition.method, IncludeOrigin = IncludeOrigin)
    all.neighbour.i <- apply(all.neighbour.i, 1, list)
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
    assign("addon.code.trace", append(get("addon.code.trace", envir = .GlobalEnv), list(sel.new.neighbour)), envir = .GlobalEnv)
    assign("addon.bic.trace", append(get("addon.bic.trace", envir = .GlobalEnv), min(neighbour.bic)), envir = .GlobalEnv)
    cat("number of add-on iterations:", as.integer(length(addon.bic.trace)), "current min BIC:", min(neighbour.bic), "\n" )
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

AddOnSearch.rxGlm <- function(merge.code,
                              varia.list,
                              model,
                              transition.method,
                              IncludeOrigin){
  current.bic <- fc.model.refit(varia.list = varia.list,
                                merge.list = merge.code,
                                mod=model)[[1]]
  num <- length(varia.list)
  all.neighbour.array <- NULL
  for (i in seq_len(num)){
    code.i <- merge.code[[i]]
    all.neighbour.i <- partition.all.neighbour(code.i, method=transition.method, IncludeOrigin = IncludeOrigin)
    all.neighbour.i <- list(unlist(apply(all.neighbour.i, 1, list), recursive = F))
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
    assign("addon.code.trace", append(get("addon.code.trace", envir = .GlobalEnv), list(sel.new.neighbour)), envir = .GlobalEnv)
    assign("addon.bic.trace", append(get("addon.bic.trace", envir = .GlobalEnv), min(neighbour.bic)), envir = .GlobalEnv)
    cat("number of add-on iterations:", as.integer(length(addon.bic.trace)), "current min BIC:", min(neighbour.bic), "\n" )
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
