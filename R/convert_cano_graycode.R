#' Convert a given graycode into its canonical format
#'
#' Convert a given graycode into its canonical format: beginning with 1, digits representing different groups are of increasing order
#'
#' @param code the input graycode which can be of any random  integer strings
#' @return the canonical graycode
#' @examples
#' a <- convert.cano.graycode(c(2, 1, 3, 4, 2)); print(a)
#' b <- sample(c(1:9), 10000, replace = T)
#' system.time({b <- convert.cano.graycode(b)})
#' print(b)
#' @export convert.cano.graycode
#'
convert.cano.graycode <- function(code){

  occur <- NULL
  uniq <- unique(code)
  update.code <- code
  for (j in (1:length(code))){
    if (code[j] %in% occur == F){
      update.code[which(code ==code[j])]<- (length(occur)+1)
      occur <- c(occur, code[j])
    }
  }
  return(update.code)
}

