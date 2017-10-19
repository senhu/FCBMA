#' Generate neighbours of a graycode with or withour pre-specified restrictions
#'
#' This algorithm randomly pick two groups and put them together, or randomly pick one group and split it into 2 groups, only generate one random neighbour each time
#'
#' @param x a numerical vector representing a graycode
#' @param method which method used to generate neighbours; either \code{"ChangeOne"} or \code{"GroupSplit"}
#' @param IncludeOrigin logical; indicating whether the input graycode should be included
#' @param group a vector indicating which elements should be grouped together
#' @param apart a vector indicating which elements should not be grouped together
#' @return a random neighbouring partition or a matrix of all neighbouring partitions
#' @author Sen Hu
#' @references Hu, S., O'Hagan, A. and Murphy, T. B. (2007). Motor Insurance Accidental Damage Claims Modelling with Factor Collapsing and Bayesian model Averaging
#' @examples
#' partition.random.neighbour(x = c(1,2,3,1), IncludeOrigin = F) # "ChangeOne" method by default
#' partition.all.neighbour(x = c(1,2,3,1), IncludeOrigin = F)
#'
#' partition.random.neighbour(c(1,2,3,4,2,4,1), method="GroupSplit", IncludeOrigin = F)
#' partition.all.neighbour(c(1,2,2,3,3,3,3), method="GroupSplit", IncludeOrigin = F)
#'
#' # the following example shows the distribution of generating random partition is not uniform
#' x <- c(1,2,2,2,3,3,3)
#' allneighbour <- NULL
#' for (j in c(1:1000)){
#'   oneneighbour <- partition.random.neighbour.GroupSplit(x)
#'   allneighbour <- rbind(allneighbour, oneneighbour)
#' }
#' dim(allneighbour)
#' uniquecode <- unique(allneighbour)
#' dim(uniquecode)
#' dim(partition.all.neighbour(x, method="GroupSplit"))
#' count1<-count2<-count3<-count4<-count5<-count6<-count7<-count8<-count9<-0
#' for (i in c(1: dim(allneighbour)[1])){
#'   if (all(uniquecode[1,] == allneighbour[i,])) count1 <- count1 + 1
#'   if (all(uniquecode[2,] == allneighbour[i,])) count2 <- count2 + 1
#'   if (all(uniquecode[3,] == allneighbour[i,])) count3 <- count3 + 1
#'   if (all(uniquecode[4,] == allneighbour[i,])) count4 <- count4 + 1
#'   if (all(uniquecode[5,] == allneighbour[i,])) count5 <- count5 + 1
#'   if (all(uniquecode[6,] == allneighbour[i,])) count6 <- count6 + 1
#'   if (all(uniquecode[7,] == allneighbour[i,])) count7 <- count7 + 1
#'   if (all(uniquecode[8,] == allneighbour[i,])) count8 <- count8 + 1
#'   if (all(uniquecode[9,] == allneighbour[i,])) count9 <- count9 + 1
#' }
#' c(count1, count2, count3, count4, count5, count6, count7, count8, count9)
#'
#' partition.all.neighbour.restricted(x = c(1,2,3,1),
#'                                    group=NULL,
#'                                    apart=c(2,3),
#'                                    IncludeOrigin = T)
#' partition.all.neighbour.restricted(x=c(1,2,3,1),
#'                                    group = c(2,3),
#'                                    apart = NULL,
#'                                    IncludeOrigin = T)
#' partition.all.neighbour.restricted(x=c(1,2,3,1),
#'                                    IncludeOrigin = T)
#'
#' @export partition.random.neighbour
#' @export partition.all.neighbour
#' @export partition.random.neighbour.ChangeOne
#' @export partition.all.neighbour.ChangeOne
#' @export partition.random.neighbour.GroupSplit
#' @export partition.all.neighbour.GroupSplit
#' @export partition.all.neighbour.restricted
#' @name partition.neighbour

NULL

#' @rdname partition.neighbour
partition.random.neighbour <- function(x, method = "ChangeOne",IncludeOrigin = FALSE){
  if (method == "ChangeOne"){return(partition.random.neighbour.ChangeOne(x, IncludeOrigin = FALSE))}
  if (method == "GroupSplit"){return(partition.random.neighbour.GroupSplit(x, IncludeOrigin = FALSE))}
}

#' @rdname partition.neighbour
partition.all.neighbour <- function(x, method = "ChangeOne",IncludeOrigin = FALSE){
  if (method == "ChangeOne"){return(partition.all.neighbour.ChangeOne(x, IncludeOrigin = FALSE))}
  if (method == "GroupSplit"){return(partition.all.neighbour.GroupSplit(x, IncludeOrigin = FALSE))}
}

#' @rdname partition.neighbour
partition.all.neighbour.restricted <- function(x, group=NULL, apart=NULL, method="ChangeOne", IncludeOrigin=FALSE){
  if (method == "ChangeOne"){
    all.neighbour <- partition.all.neighbour.ChangeOne(x, IncludeOrigin=IncludeOrigin)
  }
  if (method == "GroupSplit"){
    all.neighbour <- partition.all.neighbour.GroupSplit(x, IncludeOrigin=IncludeOrigin)
  }
  if ( is.null(group)==TRUE && is.null(apart)==TRUE ){ warning("No grouping or spliting conditions given") }
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
  new.all.neighbour <- all.neighbour[apply(all.neighbour, 1 , selection.fun),]
  return(new.all.neighbour)
}


partition.random.neighbour.ChangeOne <- function(x, IncludeOrigin = FALSE){
  which.digit <- sample(1:length(x), size=1)
  max.val <- min(  (max(as.numeric(x))+1), length(x)  )
  if (IncludeOrigin == FALSE){
    new.digit <- sample( c(1:max.val)[c(1:max.val) != x[which.digit]], size=1 )
  }
  if (IncludeOrigin == TRUE){
    new.digit <- sample( c(1:max.val), size=1 )
  }
  x[which.digit] <- new.digit
  update.x <- convert.cano.graycode(x)

  if (IncludeOrigin == FALSE){
    if (identical(update.x, x) == TRUE){
      return(Recall(x))
    }
    else {return(update.x)}
  }
  if (IncludeOrigin == TRUE){
    return(update.x)
  }
}
partition.all.neighbour.ChangeOne <- function(x, IncludeOrigin=FALSE){

  all.neighbour <- matrix(x, ncol=length(x))
  max.val <- min(  (max(as.numeric(x))+1), length(x)  )
  for (i in seq_len(length(x))){
    for (j in (seq_len(max.val)[c(1:max.val) != x[i]]) ){
      xx <- replace(x, i, j)
      update.x <- convert.cano.graycode(xx)
      if ( any(apply(all.neighbour, 1, function(y, z) isTRUE(all.equal(y, z)), update.x)) == FALSE ) {
        all.neighbour <- rbind(all.neighbour, update.x)
      }
    }
  }
  if (IncludeOrigin == FALSE){
    all.neighbour <- all.neighbour[-1,]
  }
  rownames(all.neighbour) <- NULL
  return(all.neighbour)
}
partition.random.neighbour.GroupSplit <- function(x, IncludeOrigin=FALSE){
  x <- convert.cano.graycode(x)
  if (length(unique(x)) == length(x)){step <- 1}
  if (length(unique(x)) == 1){step <- 2}
  else {step <- sample(c(1:2), 1)}

  if (step == 1){  # grouping
    if (IncludeOrigin == FALSE){ which2 <- sample(unique(x), 2, replace = F) }
    if (IncludeOrigin == TRUE){ which2 <- sample(x, 2, replace = F) }
    x[c(which(which2[1]== x), which(which2[2] == x))] <- which2[1]
    update.list.1 <- convert.cano.graycode(x)
    return(update.list.1)
  }

  if (step ==2){ # spliting
    countNum <- lapply(unique(x), function(y){sum(x==y)} )
    ToBeSelected <- unique(x)[countNum > 1]
    if (length(ToBeSelected) >= 1){
      which1 <- sample(ToBeSelected, 1)
      if (IncludeOrigin == FALSE){
        repeat{
          SplitKey <- as.factor(sample(c(1:2), sum(x == which1), replace = T))
          if (length(levels(SplitKey))>1 ) break }
        SplitGroups <- split(x[x == which1],  SplitKey)
        SplitGroups[[2]] <- rep( max(x)+1, length(SplitGroups[[2]]))
        ToBeReplaced <- unsplit(SplitGroups, SplitKey)
        ReplaceKey <- which(x == which1)
        for (i in c(1:length(ReplaceKey))){
          x[ReplaceKey[i]] <- ToBeReplaced[i]
        }
        update.list.2 <- convert.cano.graycode(x)
      }
      if (IncludeOrigin == TRUE){
        SplitKey <- as.factor(sample(c(1:2), sum(x == which1), replace = T))
        if (length(unique(SplitKey)) == 1){update.list.2 <- x}
        if (length(unique(SplitKey)) != 1){
          SplitGroups <- split(x[x == which1],  SplitKey)
          SplitGroups[[2]] <- rep( max(x)+1, length(SplitGroups[[2]]))
          ToBeReplaced <- unsplit(SplitGroups, SplitKey)
          ReplaceKey <- which(x == which1)
          for (i in c(1:length(ReplaceKey))){
            x[ReplaceKey[i]] <- ToBeReplaced[i]
          }
          update.list.2 <- convert.cano.graycode(x)
        }
      }
      return(update.list.2)
    }
    if (length(ToBeSelected) == 0){Recall(x, IncludeOrigin)}
  }
}
partition.all.neighbour.GroupSplit <- function(x, IncludeOrigin=FALSE){
  x <- convert.cano.graycode(x)
  all.neighbour <- NULL
  # grouping
  if (length(unique(x)) > 1){
    all.combn <- combinat::combn(unique(x), 2)
    combin.part <- function(whichway, y){
      y[c(which(whichway[1]== y), which(whichway[2] == y))] <- whichway[1]
      updated.list <- convert.cano.graycode(y)
      return(updated.list)
    }
    combin.neighbour <- apply(as.matrix(all.combn), 2, combin.part, y=x)
  }
  if (length(unique(x)) == 1) {combin.neighbour <- NULL}
  # spliting
  countNum <- lapply(unique(x), function(y){sum(x==y)} )
  ToBeSelected <- unique(x)[countNum > 1]
  split.part <- function(which.one, z){
    SplitKey <- as.matrix(setparts(restrictedparts(sum(z == which.one), 2))[,-1])
    split.one.part <- function(oneSplitKey, y){
      oneSplitKey <- as.factor(oneSplitKey)
      SplitGroups <- split(y[y == which.one],  oneSplitKey)
      SplitGroups[[2]] <- rep( max(y)+1, length(SplitGroups[[2]]))
      ToBeReplaced <- unsplit(SplitGroups, oneSplitKey)
      ReplaceKey <- which(y == which.one)
      for (i in c(1:length(ReplaceKey))){
        y[ReplaceKey[i]] <- ToBeReplaced[i]
      }
      updated.list <- convert.cano.graycode(y)
      return(updated.list)
    }
    updated.set <- apply(SplitKey, 2, split.one.part, y = z)
    return(updated.set)
  }
  if (length(ToBeSelected) >= 1){
    #split.neighbour <- lapply(ToBeSelected, split.part, z = x)
    split.neighbour <- NULL
    for (k in seq_len(length(ToBeSelected))){
      split.neighbour <- cbind(split.neighbour, split.part(which.one=ToBeSelected[k], z=x))
    }
  }
  if (length(ToBeSelected)==0) {split.neighbour <- NULL}
  if (IncludeOrigin == FALSE){
    return(t(as.matrix(cbind(combin.neighbour, split.neighbour))))
  }
  if (IncludeOrigin == TRUE){
    res <- t(as.matrix(cbind(combin.neighbour, split.neighbour, x)))
    rownames(res) <- NULL
    return(res)
  }
}

