#' Factor collapsing over one factor
#'
#' For a given (standard) gray code, and a variable, this function returns the collapsed version of this factor based on this given gray code
#'
#' @param variable a given variable
#' @param graycode a given graycode
#'
#' @return a collapsed version of the input variable
#'
#' @examples
#' data(sweden)
#' factor_level_collapsing(sweden$Make, c(1,2,3,4,5,6,7,8,8))
#'
#' @export

factor_level_collapsing <- function(variable,
                                    graycode){
  newvar <- variable
  level.num <- length(unique(graycode))
  for (i in c(1:level.num)){
    each_level <- c(1:length(graycode))[which(i == graycode)]
    merged_level <- levels(variable)[each_level]
    if (length(each_level) > 1){
      levels(newvar)[match(levels(variable)[each_level] , levels(newvar))] <- paste("", merged_level, sep = "", collapse=",")
    }
  }
  return(newvar)
}


