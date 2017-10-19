#' Reorganise dataset
#'
#' Given dataset, variables to be collapsed, and corresponding graycode, this function returns the collapsed dataset
#' @param data the dataset to be reorganised/collapsed
#' @param varia.list variables to be collapsed
#' @param graycode.list graycodes
#' @return new reorganised/collapsed dataset
#' @export collapse.data

collapse.data <- function(data,
                          varia.list,
                          graycode.list){
  for (i in seq_len(length(varia.list))){
    data[, which(names(data) ==varia.list[i])] <-
      factor.level.collapsing(data[, which(names(data) == varia.list[i])], graycode.list[[i]])
  }
  return(data)
}
