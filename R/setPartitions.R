#' Generate set partitions
#'
#' Enumeration of set partitions given the number of elements in a set
#'
#' @param num the number of elements in a set
#' @param group a vector indicating which elements should be grouped together
#' @param apart a vector indicating which elements should not be grouped together
#'
#' @return a matrix of all set partitions
#'
#' @note This function is based on the \code{\link[partitions]{setparts}} function from \code{"partitions"} package (R. K. S. Hankin) with modifications for FCBMA problem.
#' @references R. K. S. Hankin 2007. "Set partitions in R". Journal of Statistical Software, Volume 23, code snippet 2
#' @references R. K. S. Hankin 2006. Additive integer partitions in R. Journal of Statistical Software, Code Snippets 16(1)
#' @references Hu, S., O'Hagan, A. and Murphy, T. B. (2018). Motor insurance claims modelling with factor collapsing and Bayesian model averaging
#' @examples
#' setPartitions(4)
#' setPartitions_restricted(num=4, group=NULL, apart=c(2,3))
#' @export setPartitions
#' @export setPartitions_restricted

setPartitions <- function(num){
  Set <- as.matrix(partitions::setparts(num))
  Set <- t(apply(Set, 2, convert_canonical_graycode))
  return(Set)
}

#' @rdname setPartitions
setPartitions_restricted <- function(num, group=NULL, apart=NULL){
  Set <- setPartitions(num)
  selection.fun <- function(list){
    if ( is.null(group)==TRUE && is.null(apart)==FALSE ){
      list2 <- list[apart]
      return( length(unique(list2))==length(list2) )}
    if ( is.null(group)==FALSE && is.null(apart)==TRUE ){
      list1 <- list[group]
      return( length(unique(list1))==1 )}
    if ( is.null(group)==FALSE && is.null(apart)==FALSE ){
      list1 <- list[group]; list2 <- list[apart]
      return( length(unique(list1))==1 && length(unique(list2))==length(list2) )}
    if ( is.null(group)==TRUE && is.null(apart)==TRUE ){
      return(TRUE) }
  }
  NewSet <- Set[apply(Set, 1, selection.fun),]
  return(NewSet)
}
