#' Combining the best collapsings from factors collapsed separately
#'
#' When collapsing factors individually for the reason of computational intensity, this functions searches across the best results from each individual factor and find the best combinations across all factors
#' @param all.collapsings a list of list of all collapsing results for each factor
#' @param varia.list a vector of all factors collapsed
#' @param model full model
#' @return the best combinations
#' @examples
#' Res1 <- FactorCollapse_SA(model = mod1.1,
#'                           varia.list = c("Make"),
#'                           neighbour = "ChangeOne",
#'                           Temp = 100,
#'                           coolingRate = 0.2,
#'                           stop.Temp = 0.1,
#'                           Mtrial=20,
#'                           AddOnSearch = F)
#' Res1 <- Extract.SA.Res(Res1, threshold.val = 0.95)$State
#'
#' Res2 <- FactorCollapse_SA(model = mod1.1,
#'                           varia.list = c("Zone"),
#'                           neighbour = "ChangeOne",
#'                           Temp = 100,
#'                           coolingRate = 0.2,
#'                           stop.Temp = 0.1,
#'                           Mtrial=20,
#'                           AddOnSearch = F)
#' Res2 <- Extract.SA.Res(Res2, threshold.val = 0.95)$State
#'
#' Res3 <- FactorCollapse_SA(model = mod1.1,
#'                           varia.list = c("Kilometres"),
#'                           neighbour = "ChangeOne",
#'                           Temp = 100,
#'                           coolingRate = 0.9,
#'                           stop.Temp = 0.1,
#'                           Mtrial=30,
#'                           AddOnSearch = F)
#' Res3 <- Extract.SA.Res(Res3, threshold.val = 0.95)$State
#' z <- CombineCollapsing(all.collapsings = list(unlist(Res1$Accept.States, recursive = F),
#'                                       unlist(Res2$Accept.States, recursive = F),
#'                                       unlist(Res3$Accept.States, recursive = F)),
#'                                 varia.list=c("Make", "Zone", "Kilometres"),
#'                                 model=mod1.1)
#' ha <- combineCollapsing(all.collapsings = list(ny.1$State, ny.2$State),
#'                   varia.list = list(c("PH_GEND",
#'                                       "VEH_FUEL_TYPE",
#'                                       "VEH_TRANS",
#'                                       "DRIVEOPT",
#'                                       "CLASSOFUSE",
#'                                       "NCD_PROT",
#'                                       "STEPBACK",
#'                                       "EXCESS_TOT.factor",
#'                                       "SECONDCAR",
#'                                       "PRODUCT",
#'                                       "PAYPLAN",
#'                                       "LIC_FP_GROUP",
#'                                       "PH_PENPTS",
#'                                       "MD_LIC_CAT_A",
#'                                       "NCDLAST.factor",
#'                                       "NUMDRIVS.factor",
#'                                       "MILEAGE_B",
#'                                       "SRC_CHANNEL_A",
#'                                       "PREM_IGNITIONDISC"),
#'                                     c("COUNTYNAME") )
#'                   ,
#'                   model = train.sev.rx.1)
#' @export combineCollapsing
#'

combineCollapsing <- function(all.collapsings, varia.list, model){
  var.length <- length(all.collapsings)
  grid.states = expand.grid(all.collapsings)
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
                  bic.model.weight(sort(bic.vec, decreasing = FALSE)) )
  names(new.grid.states) <- c(varia.list, "BIC", "BMAweight")
  return(new.grid.states)
}

