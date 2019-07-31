
rxAICtoBIC <- function(rxmodel, rxdata){
  if (is.na(rxmodel$aic[1])) {stop("No AIC avaialbe for the rx-model")}
  bic.val <- rxmodel$aic - 2*rxmodel$df[1] + log(nrow(rxdata))*rxmodel$df[1]
  loglike.val <- (rxmodel$aic - 2*rxmodel$df[1])/(-2)
  return(c(bic.val, loglike.val))
}


# This tournament selection is written for function minimisation
TournamentSelection <- function(num, k, fitvec){
  pick.vec <- rep(0, num)
  for (j in c(1:num)){
    temp.min <- min(fitvec[sample(c(1:length(fitvec)), k, replace = TRUE)])
    pick.vec[j] <- which(fitvec==temp.min)[1]
  }
  return(pick.vec)
}

Gini <- function(a, p, expo = NULL, plot=FALSE) {
  if (length(a) !=  length(p)){stop("Actual and Predicted need to be equal lengths!")}
  if (!is.null(expo)) {
    temp.df <- data.frame(actual = a,
                          pred = p,
                          exposure = expo,
                          range=c(1:length(a)))
    temp.df <- temp.df[order(temp.df$pred, temp.df$range, decreasing=TRUE),]
    total.losses.act <- sum(a)
    accum.losses.act <- temp.df$actual / total.losses.act
    accum.losses.new.act <- cumsum(accum.losses.act)
    accum.losses.new.act <- c(0, accum.losses.new.act)
    total.losses.exp <- sum(expo)
    accum.losses.exp <- temp.df$exposure / total.losses.exp
    accum.losses.new.exp <- cumsum(accum.losses.exp)
    accum.losses.new.exp <- c(0, accum.losses.new.exp)
    if (plot==TRUE){
      plot(accum.losses.new.exp, accum.losses.new.act, type="l")
    }
    sum.fun <- rep(0, length(a))
    for (j in c(1:length(a))){
      sum.fun[j] <- (accum.losses.new.act[j+1]+accum.losses.new.act[j])*(accum.losses.new.exp[j+1]-accum.losses.new.exp[j])/2
    }
    gini <- sum(sum.fun)
  }
  if (is.null(expo)){
    temp.df <- data.frame(actual = a, pred = p, range=c(1:length(a)))
    temp.df <- temp.df[order(temp.df$pred, temp.df$range, decreasing=FALSE),]
    population.delta <- 1 / length(a)
    total.losses <- sum(a)
    null.losses <- rep(population.delta, length(a))
    null.losses.new <- cumsum(null.losses)
    null.losses.new <- c(0, null.losses.new)
    accum.losses <- temp.df$actual / total.losses
    accum.losses.new <- cumsum(accum.losses)
    accum.losses.new <- c(0, accum.losses.new)
    sum.fun <- rep(0, length(a))
    for (j in c(1:length(a))){
      sum.fun[j] <- (accum.losses.new[j+1]+accum.losses.new[j])*(null.losses.new[j+1]-null.losses.new[j])
    }
    gini <- 1-sum(sum.fun)
  }
  return(gini)
}
