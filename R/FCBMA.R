#' Factor collapsing (FC) with Bayesian model averaging (BMA) with lm, glm, rxGlm (package "RevoScaleR")
#'
#' This function use FC-BMA method to refit linear or generalized linear models, method chosen from complete exhaustive search, Simulated Annealing (SA), Generic Algorithm (GA).
#'
#' @param model The model used for factor collapsing.
#' @param varia.list Vector of variables to be collapsed.
#' @param method Search method to use, choose from \code{"complete"}, \code{"SA"}, \code{"GA"}.
#' @param transition.method Transition method to change one graycode to its neighbouring graycode, choose from \code{"ChangeOne"}, \code{"GroupSplit"}. The default is \code{"ChangeOne"}.
#' @param IncludeOrigin Logical; When choose neighbouring graycode, should the input graycode be included or not.
#' @param AddOnSearch Logical; used when \code{method = "SA"} or \code{"GA"}. If \code{TRUE}, after stochastic search a short greedy search will be implemented to ensure the final result is indeed global optimum. The default is \code{FALSE}.
#' @param group A (list of) vector indicating which elements must be grouped together, the default is \code{NULL}.
#' @param apart A (list of) vector indicating which elements must not be grouped together, the default is \code{NULL}.
#' @param cut Logical; if \code{FALSE}, all searched results will be presented. If \code{TRUE}, only collapsed results above a certain cumulative BMA weights threshold will be presented. The default is \code{FALSE}.
#' @param cutoff Only when \code{cut=TRUE}; threshold value of accummulated best few model weights, the default is 0.8
#' @param Temp SA parameter, initial temperature, default value is 10000.
#' @param coolingRate SA parameter, cooling rate, default value is 0.7.
#' @param stop.Temp SA parameter, stoping temperature, default value is 1e-5.
#' @param Mtrial SA parameter, number of interations at each temperature level, default is 20.
#' @param popnSize GA parameter, population size, detaulf is 20.
#' @param CrossOverRate GA parameter, crossover rate, default is 0.8.
#' @param MutationRate GA parameter, mutation rate, default is 0.1.
#' @param elitism GA parameter, elitism rate, default is 0.1.
#' @param MaxGen GA parameter, maximum generation allowed, default is 200.
#' @param verbose Logical argument controlling whether progress will be printed while the search runs. Default is \code{TRUE}.
#' @return A list with the elements
#' \describe{
#' \item{Best.Bic}{The best value of BIC.}
#' \item{Best.State}{The best partitions for all stated variables.}
#' \item{All.BIC}{All the searched BIC values.}
#' \item{All.States}{All the searched (combinations of) partitions.}
#' \item{Table}{A summary table of best few BICs, partitions, BMA weights.}
#' \item{Accept.BIC}{Accepted BIC values along SA search; only when \code{method="SA"}.}
#' \item{Accept.States}{Accepted (combinations of) partitions along SA search; only when \code{method="SA"}.}
#' \item{Best.BIC.trace}{Trace of best BIC values along GA search; only when \code{method="GA"}.}
#' \item{AddOn.BestBIC}{The best BIC value found through the add-on greedy search, if \code{AddOnSearch=TRUE}.}
#' \item{AddOn.BestState}{The best (combination of) partition found through the add-on greedy search, if \code{AddOnSearch=TRUE}.}
#' \item{AddOn.BIC.trace}{The best BIC values found in each iteration through the add-on greedy search, if \code{AddOnSearch=TRUE}.}
#' \item{AddOn.State.trace}{The best (combination of) partitions found in each iteration through the add-on greedy search, if \code{AddOnSearch=TRUE}.}
#' }
#' @note This version of stochastic search changes on the simulated annealing part. It changed the stochastic search method at each iteration, only one variable is picked and change to one of its neighbours the probability of being picked is based on the number of levels You also have the option of which method to select a neighbour of a graycode
#' @examples
#' # get the Sweden TP insurance data set
#' library(faraway)
#' data("motorins")
#' dataset <- motorins
#' dataset$Bonus <- as.factor(dataset$Bonus)
#' dataset$Kilometres <- factor(dataset$Kilometres, ordered=F)
#'
#' # Frequency Poisson GLM
#' freq1 <- glm(Claims~Zone+Bonus+Make+Kilometres,
#'              offset=log(Insured),
#'              data=dataset,
#'              family="poisson")
#' summary(freq1)
#'
#' freq.comp <- FCBMA(freq1,
#'                     varia.list = c("Kilometres"),
#'                     method = "complete",
#'                     verbose = TRUE)
#' system.time({
#'   freq.sa <- FCBMA(model=freq1,
#'                    varia.list = c("Kilometres","Make"),
#'                    method = "SA",
#'                    transition.method = "ChangeOne",
#'                    AddOnSearch = TRUE,
#'                    Temp = 1000,
#'                    coolingRate = 0.6,
#'                    stop.Temp = 1e-5,
#'                    Mtrial = 20,
#'                    verbose = TRUE)
#' })
#' freq.sa$Best.State
#' freq.sa$Best.BIC
#' system.time({
#'   freq.ga <- FCBMA(model=freq1,
#'                    varia.list = c("Make"),
#'                    method = "GA",
#'                    transition.method = "ChangeOne",
#'                    popnSize=20,CrossOverRate=0.8,MutationRate=0.1,
#'                    elitism=0.1,
#'                    MaxGen=100,
#'                    verbose = TRUE)
#' })
#' freq.ga$Best.State
#' freq.ga$Best.BIC
#'
#' @export FCBMA

FCBMA <-function(model,
                  varia.list,
                  method = c("complete", "SA", "GA"),
                  transition.method="ChangeOne",
                  IncludeOrigin = FALSE,
                  AddOnSearch = FALSE,
                  group=NULL, apart=NULL,
                  cut = FALSE,cutoff=0.8,
                  Temp = 10000,coolingRate = 0.7,stop.Temp = 1e-5,Mtrial = 20,
                  popnSize=20,CrossOverRate=0.8,MutationRate=0.1,elitism=0.1,MaxGen=200,
                  verbose=TRUE
                  ){
  if (method == "complete"){
    newdat <- eval(getCall(model)$data,envir = environment(formula(model)))
    nvar <- length(varia.list)
    level.num <- NULL
    collapse.comb <- NULL
    comb.num <- NULL
    if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
      for (i in seq_len(nvar)){
        var.temp <- newdat[,which(colnames(newdat) == varia.list[i])]
        n.temp <- nlevels(var.temp)
        level.num <- c(level.num, n.temp)
        AA <- setPartitions(n.temp)
        collapse.comb <- append(collapse.comb, list(AA))
        comb.num <- c(comb.num, nrow(AA))
      }
    } else {
      for (i in seq_len(nvar)){
        var.temp <- newdat[,which(colnames(newdat) == varia.list[i])]
        n.temp <- nlevels(var.temp)
        level.num <- c(level.num, n.temp)
        AA <- setPartitions.restricted(n.temp, group=group[[i]], apart=apart[[i]])
        collapse.comb <- append(collapse.comb, list(AA))
        comb.num <- c(comb.num, nrow(AA))
      }
    }
    search.grid <- as.matrix(expand.grid(lapply(comb.num, seq_len)) )
    names(search.grid) <- NULL
    BICval.vec <- NULL
    Graycode.vec <- NULL
    Partition.vec <- NULL
    Graycode.mat <- NULL
    for (a in c(1:nrow(search.grid))){
      merge.list <- NULL
      for (t in seq_len(ncol(search.grid))){
        merge.list <- append(merge.list, list(collapse.comb[[t]][search.grid[a, t], ]))
      }
      BICval.vec <- c(BICval.vec, fc.model.refit(varia.list = varia.list,
                                                 merge.list = merge.list,
                                                 mod=model)[[1]]  )
      Graycode.mat <- append(Graycode.mat, list(merge.list))
      r1 <- NULL; r2 <- NULL
      for (i in c(1:length(merge.list))){
        r1 <- cbind(r1, paste(c("(",paste(merge.list[[i]], collapse = ","),")"),collapse = ""))
        r2 <- cbind(r2, graycode.to.partition(merge.list[[i]]))
      }
      Graycode.vec <- rbind(Graycode.vec, r1)
      Partition.vec <- rbind(Partition.vec, r2)
      if (verbose){print(a)}
    }
    bestbic <- min(BICval.vec)
    beststate <- Graycode.mat[[which.min(BICval.vec)]]
    BMAweight.vec <- bic.model.weight(BICval.vec)
    Table<- as.data.frame(cbind(Graycode.vec, Partition.vec, as.numeric(BICval.vec), as.numeric(BMAweight.vec)))
    names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
    Table <- Table[order(Table$BIC, decreasing=F), ]
    rownames(Table) <- NULL
    if (cut == FALSE){return(list(Best.BIC=bestbic,
                                  Best.State=beststate,
                                  All.BIC=BICval.vec,
                                  All.States=Graycode.mat,
                                  Table=Table))}
    if (cut == TRUE) {
      num.cutoff <- which(cumsum(sort(BMAweight.vec, decreasing = TRUE)) >= cutoff)[1]
      return(list(Best.BIC=bestbic,
                  Best.State=beststate,
                  All.BIC=BICval.vec,
                  All.States=Graycode.mat,
                  Table[1:num.cutoff,]))}
  }
  if (method == "SA"){
    newdat <- eval(getCall(model)$data, envir = environment(formula(model)))
    nvar <- length(varia.list)
    level.num <- NULL
    for (i in seq_len(nvar)){
      var.temp <- newdat[,which(colnames(newdat) == varia.list[i])]
      n.temp <- nlevels(var.temp)
      level.num <- c(level.num, n.temp)
    }
    All_states <- NULL
    All_bic <- NULL
    Accept_states <- NULL
    Accept_bic <- NULL
    S_c_list <- NULL
    S_c_mat <- NULL
    for (j in seq_len(nvar)){
      if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
        S_c_temp <- convert.cano.graycode( sample(c(1:level.num[j]), size = level.num[j], replace=TRUE) )
      } else {svec <- sample(c(1:level.num[j]), size = level.num[j], replace=TRUE)
              svec[group[[j]]] <- svec[group[[j]]][1]
              svec[apart[[j]]] <- seq_len(length(apart[[j]]))
              S_c_temp <- convert.cano.graycode(svec)}
      S_c_list <- append(S_c_list, list(S_c_temp))
    }
    F_b <- F_n <- F_c <-  fc.model.refit(varia.list = varia.list,
                                         merge.list = S_c_list,
                                         mod = model)[[1]]
    All_states <- append(All_states, list(S_c_mat))
    All_bic <- c(All_bic, F_c)
    accepted.vec <- NULL
    if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
      while (Temp > stop.Temp){
        M <- 1
        accepted.M <- 1
        while (M <= Mtrial ){
          which_var <- sample(c(1:nvar),1, prob=level.num)
          first.all.neighbour <- partition.all.neighbour(x = S_c_list[[which_var]], method=transition.method, IncludeOrigin = IncludeOrigin)
          first.neighbour.num <- sample(c(1:nrow(first.all.neighbour)), 1, replace=T)
          S_n_list <- S_c_list
          S_n_list[[which_var]] <- first.all.neighbour[first.neighbour.num, ]
          F_n <- fc.model.refit(varia.list=varia.list,
                                merge.list=S_n_list,
                                mod = model)[[1]]
          All_states <- append(All_states, list(S_n_list))
          All_bic <- c(All_bic, F_n)
          if (F_n < F_c || runif(1, 0, 1) < exp( (F_c - F_n) / Temp)) {
            Accept_states <- append(Accept_states, list(S_n_list))
            Accept_bic <- c(Accept_bic, F_n)
            F_c <- F_n
            S_c_list <- S_n_list
            accepted.M <- accepted.M + 1
          }
          if (F_n < F_b) {
            S_best_list <- S_n_list
            F_b <- F_n
          }
          M <- M+1
          if (verbose){print(c(M, Temp))}
        }
        accepted.vec <- c(accepted.vec, accepted.M)
        Temp <- Temp * coolingRate
        #Temp <- Temp * (1 + (F_c - F_b)/F_c )
      }
    } else {
      while (Temp > stop.Temp){
        M <- 1
        accepted.M <- 1
        while (M <= Mtrial ){
          which_var <- sample(c(1:nvar),1, prob=level.num)
          first.all.neighbour <- partition.all.neighbour.restricted(x = S_c_list[[which_var]],
                                                                    group=group[[which_var]],
                                                                    apart=apart[[which_var]],
                                                                    method=transition.method,
                                                                    IncludeOrigin = IncludeOrigin)
          first.neighbour.num <- sample(c(1:nrow(first.all.neighbour)), 1, replace=T)
          S_n_list <- S_c_list
          S_n_list[[which_var]] <- first.all.neighbour[first.neighbour.num, ]
          F_n <- fc.model.refit(varia.list=varia.list,
                                merge.list=S_n_list,
                                mod = model)[[1]]
          All_states <- append(All_states, list(S_n_list))
          All_bic <- c(All_bic, F_n)
          if (F_n < F_c || runif(1, 0, 1) < exp( (F_c - F_n) / Temp)) {
            Accept_states <- append(Accept_states, list(S_n_list))
            Accept_bic <- c(Accept_bic, F_n)
            F_c <- F_n
            S_c_list <- S_n_list
            accepted.M <- accepted.M + 1
          }
          if (F_n < F_b) {
            S_best_list <- S_n_list
            F_b <- F_n
          }
          M <- M+1
          if (verbose){print(c(M, Temp))}
        }
        accepted.vec <- c(accepted.vec, accepted.M)
        Temp <- Temp * coolingRate
        #Temp <- Temp * (1 + (F_c - F_b)/F_c )
      }
    }
    if (cut){
      take.n <- which(cumsum(bic.model.weight( sort(unique(Accept_bic), decreasing=F) )) >= cutoff)[1]
    } else {take.n <- length(unique(Accept_bic))}
    take.bic <- sort(unique(Accept_bic))[1:take.n]
    take.weight <- bic.model.weight( sort(unique(Accept_bic), decreasing=F) )[1:take.n]
    take.state <- NULL
    for (i in seq_len(take.n)){
      temp.state <- unlist(Accept_states[ which( Accept_bic == sort(unique(Accept_bic))[i] )[1] ], recursive = F)
      take.state <- append(take.state, list(temp.state))
    }
    Table <- list(N=take.n,
                  BIC=take.bic,
                  Weights=take.weight,
                  State=take.state)
    if (AddOnSearch){
      print("Start of add-on greedy search")
      assign("addon.bic.trace", NULL, envir = .GlobalEnv)
      assign("addon.code.trace", NULL, envir = .GlobalEnv)
      res.add <- AddOnSearch(merge.code = S_best_list,
                             varia.list = varia.list,
                             model = model,
                             transition.method=transition.method,
                             IncludeOrigin = IncludeOrigin)
      AddOn.BestBIC <- res.add[[1]]
      AddOn.BestState <- res.add[[2]]
      AddOn.BIC.trace <- get("addon.bic.trace", envir = .GlobalEnv)
      AddOn.State.trace <- get("addon.code.trace", envir = .GlobalEnv)

      return(list(Best.BIC = F_b,
                  Best.State = S_best_list,
                  All.BIC = All_bic,
                  All.States = All_states,
                  Accept.BIC = Accept_bic,
                  Accept.States = Accept_states,
                  Table = Table,
                  AddOn.BestBIC = AddOn.BestBIC,
                  AddOn.BestState = AddOn.BestState,
                  AddOn.BIC.trace = AddOn.BIC.trace,
                  AddOn.State.trace = AddOn.State.trace
                  #accepted.vec = accepted.vec
      ))
    } else{
      return(list(Best.BIC = F_b,
                  Best.State = S_best_list,
                  All.BIC = All_bic,
                  All.States = All_states,
                  Accept.BIC = Accept_bic,
                  Accept.States = Accept_states,
                  Table = Table
                  #accepted.vec = accepted.vec
      ))
    }
  }
  if (method == "GA"){
    if (AddOnSearch){
      assign("addon.bic.trace", NULL, envir = .GlobalEnv)
      assign("addon.code.trace", NULL, envir = .GlobalEnv)
    }
    newdat <- eval(getCall(model)$data, envir = environment(formula(model)))
    nvar <- length(varia.list)
    level.num <- NULL
    for (i in seq_len(nvar)){
      var.temp <- newdat[,which(colnames(newdat) == varia.list[i])]
      n.temp <- nlevels(var.temp)
      level.num <- c(level.num, n.temp)
    }
    generation.k <- 0
    Popn <- NULL
    if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
      for (each in c(1:popnSize)){
        S_c_list <- NULL
        for (j in seq_len(nvar)){
          S_c_temp <- convert.cano.graycode( sample(c(1:level.num[j]), size = level.num[j], replace=TRUE) )
          S_c_list <- append(S_c_list, list(S_c_temp))
        }
        Popn <- append(Popn, list(S_c_list))
      }
    } else {
      for (each in c(1:popnSize)){
        S_c_list <- NULL
        for (j in seq_len(nvar)){
          svec <- sample(c(1:level.num[j]), size = level.num[j], replace=TRUE)
          svec[group[[j]]] <- svec[group[[j]]][1]
          svec[apart[[j]]] <- seq_len(length(apart[[j]]))
          S_c_temp <- convert.cano.graycode(svec)
          S_c_list <- append(S_c_list, list(S_c_temp))
        }
        Popn <- append(Popn, list(S_c_list))
      }
    }
    fitness.calc.function <- function(onelist){
      temp.res <- fc.model.refit(varia.list= varia.list ,
                                 merge.list= onelist,
                                 mod=model)[[1]]
      return(temp.res)
    }
    fitness.vec <- sapply(Popn, fitness.calc.function)
    fittestOne <- min(fitness.vec)

    # internal functions
    Crossover.function <- function(listmat){
      if ( length(listmat) %% 2 !=0 ){stop("The number of parents is not even number!")}
      if ( length(listmat) %% 2 == 0){
        offspring.mat <- NULL
        if (nvar==1){
          for ( i in seq(1, length(listmat), by=2) ){
            parent.1 <- unlist(listmat[i])
            parent.2 <- unlist(listmat[i+1])
            cross.point.1 <- sort( sample( c(1:(level.num-1)), 1, replace=FALSE), decreasing=FALSE)
            mask.1 <- c(rep(TRUE, cross.point.1), rep(FALSE, sum(level.num)-cross.point.1))
            mask.2 <- c(rep(FALSE, cross.point.1), rep(TRUE, sum(level.num)-cross.point.1))
            #mask.3 <- c(rep(FALSE, cross.point.1), rep(FALSE, (cross.point.2-cross.point.1)), rep(TRUE, sum(level.num)-cross.point.2))
            offspring.1 <- list(c(parent.1[mask.1], parent.2[mask.2]))
            offspring.2 <- list(c(parent.2[mask.1], parent.1[mask.2]))
            offspring.mat <- append(offspring.mat, list(offspring.1))
            offspring.mat <- append(offspring.mat, list(offspring.2))
          }
        } else {
          for ( i in seq(1, length(listmat), by=2) ){
            parent.1 <- unlist(listmat[i])
            parent.2 <- unlist(listmat[i+1])
            cross.point <- sort( sample( c(1:(nvar-1)), 1, replace=FALSE), decreasing=FALSE)
            cross.point.1 <- cumsum(level.num)[cross.point[1]]
            #cross.point.2 <- cumsum(level.num)[cross.point[2]]
            mask.1 <- c(rep(TRUE, cross.point.1), rep(FALSE, sum(level.num)-cross.point.1))
            mask.2 <- c(rep(FALSE, cross.point.1), rep(TRUE, sum(level.num)-cross.point.1))
            #mask.3 <- c(rep(FALSE, cross.point.1), rep(FALSE, (cross.point.2-cross.point.1)), rep(TRUE, sum(level.num)-cross.point.2))
            offspring.1 <- c(parent.1[mask.1], parent.2[mask.2])
            offspring.2 <- c(parent.2[mask.1], parent.1[mask.2])
            offspring.1.list <- SplitAt(offspring.1, cumsum(level.num))
            offspring.2.list <- SplitAt(offspring.2, cumsum(level.num))
            offspring.mat <- append(offspring.mat, list(offspring.1.list))
            offspring.mat <- append(offspring.mat, list(offspring.2.list))
          }
        }
        return(offspring.mat)
      }
    }
    Mutation.function <- function(thelist){
      which.var <- sample( c(1:nvar),1, prob=level.num )
      collapse.old <- unlist(thelist[[which.var]])
      all.neighbour <- partition.all.neighbour(x = collapse.old, method = transition.method, IncludeOrigin = IncludeOrigin)
      collapse.new <- all.neighbour[ sample(c(1:nrow(all.neighbour)), 1) , ]
      thelist[[which.var]] <- collapse.new
      return(thelist)
    }
    Mutation.function.restricted <- function(thelist){
      which.var <- sample( c(1:nvar),1, prob=level.num )
      collapse.old <- unlist(thelist[[which.var]])
      all.neighbour <- partition.all.neighbour.restricted(x = collapse.old,
                                                          group=group[[which.var]],
                                                          apart=apart[[which.var]],
                                                          method = transition.method,
                                                          IncludeOrigin = IncludeOrigin)
      collapse.new <- all.neighbour[ sample(c(1:nrow(all.neighbour)), 1) , ]
      thelist[[which.var]] <- collapse.new
      return(thelist)
    }
    if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
      while (generation.k <= MaxGen){
        # elitism
        elitism.num <- round(popnSize*elitism)
        NewPopn <- Popn[order(fitness.vec, decreasing=FALSE)[1:elitism.num]]
        Popn[ order(fitness.vec, decreasing=FALSE)[1:elitism.num] ] <-NULL
        OldPopn <- Popn
        OldFitness <- fitness.vec[-(order(fitness.vec, decreasing=FALSE)[1:elitism.num])]
        # copy # using tounament rule
        copy.popn <- OldPopn[TournamentSelection(round((1- CrossOverRate)*length(OldPopn)), k=2, OldFitness) ]
        # using roulette rule #copy.popn <- OldPopn[sample(c(1:length(OldFitness)), size = round((1- CrossOverRate)*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
        NewPopn <- append(NewPopn, copy.popn)
        # cross-over # one point cross over first using tournament rule
        cross.mat <- OldPopn[TournamentSelection(round(CrossOverRate*length(OldPopn)), k=2, OldFitness)]
        # using roulette rule
        #cross.mat <- OldPopn[sample(c(1:length(OldFitness)), size = round(CrossOverRate*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
        crossover.popn <- Crossover.function(cross.mat)
        NewPopn <- append(NewPopn, crossover.popn)
        # mutation
        mutation.vec <- sample(c( (round(elitism.num)+1):(length(NewPopn)) ), round(MutationRate*popnSize), replace = FALSE)
        mutation.mat <- NewPopn[mutation.vec]
        NewPopn[mutation.vec] <- NULL
        mutation.popn <- lapply(mutation.mat, Mutation.function)
        NewPopn <- append(NewPopn, mutation.popn)
        fitness.vec <- sapply(NewPopn, fitness.calc.function)
        fittestOne <- c(fittestOne, min(fitness.vec))
        Popn <- NewPopn
        generation.k <- generation.k + 1
        if (verbose) {print(generation.k)}
      }
    } else {
      while (generation.k <= MaxGen){
        # elitism
        elitism.num <- round(popnSize*elitism)
        NewPopn <- Popn[order(fitness.vec, decreasing=FALSE)[1:elitism.num]]
        Popn[ order(fitness.vec, decreasing=FALSE)[1:elitism.num] ] <-NULL
        OldPopn <- Popn
        OldFitness <- fitness.vec[-(order(fitness.vec, decreasing=FALSE)[1:elitism.num])]
        # copy # using tounament rule
        copy.popn <- OldPopn[TournamentSelection(round((1- CrossOverRate)*length(OldPopn)), k=2, OldFitness) ]
        # using roulette rule #copy.popn <- OldPopn[sample(c(1:length(OldFitness)), size = round((1- CrossOverRate)*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
        NewPopn <- append(NewPopn, copy.popn)
        # cross-over # one point cross over first using tournament rule
        cross.mat <- OldPopn[TournamentSelection(round(CrossOverRate*length(OldPopn)), k=2, OldFitness)]
        # using roulette rule
        #cross.mat <- OldPopn[sample(c(1:length(OldFitness)), size = round(CrossOverRate*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
        crossover.popn <- Crossover.function(cross.mat)
        NewPopn <- append(NewPopn, crossover.popn)
        # mutation
        mutation.vec <- sample(c( (round(elitism.num)+1):(length(NewPopn)) ), round(MutationRate*popnSize), replace = FALSE)
        mutation.mat <- NewPopn[mutation.vec]
        NewPopn[mutation.vec] <- NULL
        mutation.popn <- lapply(mutation.mat, Mutation.function.restricted)
        NewPopn <- append(NewPopn, mutation.popn)
        fitness.vec <- sapply(NewPopn, fitness.calc.function)
        fittestOne <- c(fittestOne, min(fitness.vec))
        Popn <- NewPopn
        generation.k <- generation.k + 1
        if (verbose) {print(generation.k)}
      }
    }
    Best.State <- unlist(Popn[which.min(fitness.vec)], recursive = FALSE)
    Best.BIC <- min(fitness.vec)
    All.Popn.States <- Popn
    All.BIC <- fitness.vec
    Best.BIC.trace <- fittestOne
    if (cut){ take.n <- which(cumsum(bic.model.weight( sort(unique(Accept_bic), decreasing=F) )) >= cutoff)[1] }
    else {take.n <- length(unique(Accept_bic))}
    take.bic <- sort(unique(All.BIC))[1:take.n]
    take.weight <- bic.model.weight( sort(unique(All.BIC), decreasing=F) )[1:take.n]
    take.state <- NULL
    for (i in seq_len(take.n)){
      temp.state <- unlist(All.Popn.States[ which(All.BIC==sort(unique(All.BIC))[i] )[1] ], recursive = F)
      take.state <- append(take.state, list(temp.state))
    }
    Table <- list(N=take.n,
                  BIC=take.bic,
                  Weights=take.weight,
                  State=take.state)
    if (AddOnSearch==TRUE){
      print("Start of add-on greedy search")
      assign("addon.bic.trace", NULL, envir = .GlobalEnv)
      assign("addon.code.trace", NULL, envir = .GlobalEnv)
      res.add <- AddOnSearch(merge.code = Best.State,
                             varia.list = varia.list,
                             model = model,
                             transition.method=transition.method,
                             IncludeOrigin = IncludeOrigin)
      AddOn.BestBIC <- res.add[[1]]
      AddOn.BestState <- res.add[[2]]
      AddOn.BIC.trace <- get("addon.bic.trace", envir = .GlobalEnv)
      AddOn.State.trace <- get("addon.code.trace", envir = .GlobalEnv)
      return(list(Best.BIC = Best.BIC,
                  Best.State = Best.State,
                  All.BIC = All.BIC,
                  All.Popn.States = All.Popn.States,
                  Best.BIC.trace = Best.BIC.trace,
                  AddOn.BestBIC = AddOn.BestBIC,
                  AddOn.BestState = AddOn.BestState,
                  AddOn.BIC.trace = AddOn.BIC.trace,
                  AddOn.State.trace = AddOn.State.trace))
    } else{
      return(list(Best.BIC = Best.BIC,
                  Best.State = Best.State,
                  All.BIC = All.BIC,
                  All.Popn.States = All.Popn.States,
                  Best.BIC.trace = Best.BIC.trace ))
    }
  }
  else {stop('Unknown method. Please use "complete", "SA" or "GA" ')}

}
