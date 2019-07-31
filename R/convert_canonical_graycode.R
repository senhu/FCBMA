#' Convert a given graycode into its canonical format
#'
#' Convert a given graycode into its canonical format: beginning with 1, digits representing different groups are of increasing order
#'
#' @param code An interger vector, i.e. the input graycode, which can be of any random integer vectors
#' @return An integer vector of the canonical graycode
#' @examples
#' convert_canonical_graycode(c(2, 1, 3, 4, 2))
#' a <- sample(c(1:9), 1000, replace = TRUE)
#' convert_canonical_graycode(a)
#'
#' @export

convert_canonical_graycode <- function(code){
  occur <- NULL
  uniq <- unique(code)
  update.code <- code
  for (j in (1:length(code))){
    if (code[j] %in% occur == F){
      update.code[which(code ==code[j])] <- (length(occur)+1)
      occur <- c(occur, code[j])
    }
  }
  return(update.code)
}

