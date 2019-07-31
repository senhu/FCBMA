#' Combining the best collapsings from factors collapsed separately
#'
#' When collapsing factors individually for the reason of computational intensity, this functions searches across the best results from each individual factor and find the best combinations across all factors
#'
#' @param all.collapsings a list of list of all collapsing results for each factor
#' @param varia.list a vector of all factors collapsed
#' @param model full model
#'
#' @return the best combinations
#'
#' @examples
#' \donttest{
#' r1 <- FCBMA(model = mod1,
#'             varia.list = c("Make"),
#'             method = "SA",
#'             transition.method = "ChangeOne")
#' res1 <- Extract.SA.Res(r1, threshold.val = 0.95)$State
#'
#' r2 <- FCBMA(model = mod1,
#'             varia.list = c("Zone"),
#'             method = "SA",
#'             transition.method = "ChangeOne")
#' res2 <- Extract.SA.Res(r2, threshold.val = 0.95)$State
#'
#' r3 <- FCBMA(model = mod1,
#'             varia.list = c("Kilometres"),
#'             method = "SA",
#'             transition.method = "ChangeOne")
#' res3 <- Extract.SA.Res(r3, threshold.val = 0.95)$State
#' z <- CombineCollapsing(all.collapsings = list(unlist(Res1$Accept.States, recursive = F),
#'                                               unlist(Res2$Accept.States, recursive = F),
#'                                               unlist(Res3$Accept.States, recursive = F)),
#'                                 varia.list=c("Make", "Zone", "Kilometres"),
#'                                 model=mod1)
#' }
#'
#' @export combineCollapsing

combineCollapsing <- function(all.collapsings,
                              varia.list,
                              model){
  var.length <- length(all.collapsings)
  grid.states <- expand.grid(all.collapsings)
  temp.func <- function(row){
    temp.list <- unlist(unlist(as.list(row), recursive = F), recursive = F)
    temp.bic <- fc.model.refit(varia.list = unlist(varia.list),
                               merge.list = temp.list,
                               mod = model)[[1]]
    return(temp.bic)
  }
  bic.vec <- NULL
  for (i in seq_len(dim(grid.states)[1])){bic.vec <- append(bic.vec, temp.func(grid.states[i,]))}
  new.grid.states <- cbind(grid.states[order(bic.vec, decreasing = FALSE),],
                  sort(bic.vec, decreasing = FALSE),
                  BMAweight_bic(sort(bic.vec, decreasing = FALSE)) )
  names(new.grid.states) <- c(varia.list, "BIC", "BMAweight")
  return(new.grid.states)
}

