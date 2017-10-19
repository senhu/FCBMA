#' Tournament Selection
#'
#' Tournament Selection
#' @param num number of values required
#' @param k number of selection each time
#' @param fitvec vector of fitness value of all elements to be selected
#' @return positions of elements selected
#' @note This tournament selection is written for function minimisation
#' @examples
#' fitvec <- rnorm(10000, 0, 100)
#' TournamentSelection(200, 5, fitvec)
#' @export TournamentSelection

TournamentSelection <- function(num, k, fitvec){
  pick.vec <- rep(0, num)
  for (j in c(1:num)){
    temp.min <- min(fitvec[sample(c(1:length(fitvec)), k, replace = TRUE)])
    pick.vec[j] <- which(fitvec==temp.min)[1]
  }
  return(pick.vec)
}

