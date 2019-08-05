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
#' @param control.SA A list of control parameters for SA.
#' \describe{
#' \item{init.Temp}{Initial temperature, default value is 10000.}
#' \item{stop.Temp}{Stoping temperature, default value is 1e-5.}
#' \item{coolingRate}{Cooling rate, default value is 0.7.}
#' \item{Mtrial}{Number of interations at each temperature level, default is 20.}
#' }
#' @param control.GA A list of control parameters for GA.
#' \describe{
#' \item{popnSize}{Population size, detaulf is 20.}
#' \item{CrossOverRate}{Crossover rate, default is 0.8.}
#' \item{MutationRate}{Mutation rate, default is 0.1.}
#' \item{elitism}{Elitism rate, default is 0.1.}
#' \item{MaxGen}{Maximum generation allowed, default is 200.}
#' }
#' @param verbose Logical argument controlling whether progress will be printed while the search runs. Default is \code{TRUE}.
#'
#' @return A list with the elements
#' \describe{
#' \item{method}{The search method used in finding the optimal factor collapsings.}
#' \item{variable}{The variables collapsed.}
#' \item{best.bic}{The best value of BIC, corresponding to the \code{best.state}.}
#' \item{best.state}{The best partitions of the stated variables, in graycode format.}
#' \item{all.bic}{All the searched BIC values, in an ascending order.}
#' \item{all.states}{A list, all the searched (combinations of) partitions, in graycode format, corresponding to \code{all.bic}.}
#' \item{all.partitions}{A list, sll the searched (combinations of) partitions, corresponding to \code{all.bic}.}
#' \item{all.weights}{All the BMA weights of searched partitions, , corresponding to \code{all.bic}, in a descending order.}
#' \item{table}{A summary table of best few BICs, partitions, BMA weights.}
#' \item{model}{The input baseline model based on which factor collapsing is implemented.}
#' \item{addon}{An indicator of whether add-on greedy search is utilised.}
#' }
#'
#' @examples
#' data("sweden")
#' m1 <- glm(Claims ~ Kilometres+Zone+Bonus+Make, offset = log(Insured),
#'           data = sweden, family = "poisson")
#' summary(m1)
#'
#' # complete search
#' m2 <- FCBMA(model = m1,
#'             varia.list = c("Kilometres"),
#'             method = "complete",
#'             verbose = FALSE)
#'
#' \donttest{
#' # SA search
#' m3 <- FCBMA(model = m1,
#'             varia.list = c("Kilometres","Make"),
#'             method = "SA",
#'             transition.method = "ChangeOne",
#'             AddOnSearch = TRUE,
#'             control.SA = list(init.Temp = 1000,
#'                               coolingRate = 0.6,
#'                               stop.Temp = 1e-5,
#'                               Mtrial = 20),
#'             verbose = TRUE)
#' m3$best.state
#' m3$best.bic
#'
#' # GA search
#' m4 <- FCBMA(model=m1,
#'             varia.list = c("Kilometres", "Make"),
#'             method = "GA",
#'             transition.method = "ChangeOne",
#'             AddOnSearch= TRUE,
#'             control.GA = list(popnSize=20,
#'                               CrossOverRate=0.8,
#'                               MutationRate=0.1,
#'                               elitism=0.1,
#'                               MaxGen=100),
#'             verbose = TRUE)
#' m4$best.state
#' m4$best.bic
#' }
#'
#' @export FCBMA

FCBMA <- function(model,
                  varia.list,
                  method = "complete", #c("complete", "SA", "GA"),
                  transition.method="ChangeOne",
                  IncludeOrigin = FALSE,
                  AddOnSearch = FALSE,
                  group=NULL, apart=NULL,
                  # SA parameters
                  control.SA = list(init.Temp = 10000,
                                    stop.Temp = 1e-5,
                                    coolingRate = 0.7,
                                    Mtrial = 20),
                  # GA parameters
                  control.GA = list(popnSize = 20,
                                    CrossOverRate = 0.8,
                                    MutationRate = 0.1,
                                    elitism = 0.1,
                                    MaxGen = 200),
                  verbose=TRUE
                  ){
  switch(method,
         complete = {
           newdat <- eval(stats::getCall(model)$data,envir = environment(stats::formula(model)))
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
               AA <- setPartitions_restricted(n.temp, group=group[[i]], apart=apart[[i]])
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
               r1 <- cbind(r1, merge.list[i])
               # r1 <- cbind(r1, paste(c("(",paste(merge.list[[i]], collapse = ","),")"),collapse = ""))
               r2 <- cbind(r2, graycode_to_partition(merge.list[[i]]))
             }
             Graycode.vec <- rbind(Graycode.vec, r1)
             Partition.vec <- rbind(Partition.vec, r2)
             if (verbose){print(a)}
           }
           bestbic <- min(BICval.vec)
           beststate <- Graycode.mat[[which.min(BICval.vec)]]
           BMAweight.vec <- BMAweight_bic(BICval.vec)
           Table<- as.data.frame(cbind(Graycode.vec,
                                       Partition.vec,
                                       as.numeric(BICval.vec),
                                       as.numeric(BMAweight.vec)))
           names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
           indx <- order(BICval.vec, decreasing=FALSE)
           Table <- Table[indx, ]
           rownames(Table) <- NULL
           result <- list(method = "complete",
                          variable = varia.list,
                          best.state = beststate,
                          best.bic = bestbic,
                          all.bic = BICval.vec[indx],
                          all.states = Graycode.vec[indx],
                          all.weights = BMAweight.vec[indx],
                          all.partitions = Partition.vec[indx],
                          table = Table,
                          model = model)
           structure(result, class = c('FCBMA'))
         },
         SA = {
           newdat <- eval(stats::getCall(model)$data, envir = environment(stats::formula(model)))
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
               S_c_temp <- convert_canonical_graycode( sample(c(1:level.num[j]), size = level.num[j], replace=TRUE) )
             } else {svec <- sample(c(1:level.num[j]), size = level.num[j], replace=TRUE)
             svec[group[[j]]] <- svec[group[[j]]][1]
             svec[apart[[j]]] <- seq_len(length(apart[[j]]))
             S_c_temp <- convert_canonical_graycode(svec)}
             S_c_list <- append(S_c_list, list(S_c_temp))
           }
           F_b <- F_n <- F_c <-  fc.model.refit(varia.list = varia.list,
                                                merge.list = S_c_list,
                                                mod = model)[[1]]
           All_states <- append(All_states, list(S_c_mat))
           All_bic <- c(All_bic, F_c)
           #accepted.vec <- NULL
           if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
             Temp <- control.SA$init.Temp
             while (Temp > control.SA$stop.Temp){
               M <- 1
               accepted.M <- 1
               while (M <= control.SA$Mtrial ){
                 which_var <- sample(c(1:nvar),1, prob=level.num)
                 first.all.neighbour <- partition_all_neighbour(x = S_c_list[[which_var]],
                                                                method=transition.method,
                                                                IncludeOrigin = IncludeOrigin)
                 first.neighbour.num <- sample(c(1:nrow(first.all.neighbour)), 1, replace=TRUE)
                 S_n_list <- S_c_list
                 S_n_list[[which_var]] <- first.all.neighbour[first.neighbour.num, ]
                 F_n <- fc.model.refit(varia.list=varia.list,
                                       merge.list=S_n_list,
                                       mod = model)[[1]]
                 All_states <- append(All_states, list(S_n_list))
                 All_bic <- c(All_bic, F_n)
                 if (F_n < F_c || stats::runif(1, 0, 1) < exp( (F_c - F_n) / Temp)) {
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
               #accepted.vec <- c(accepted.vec, accepted.M)
               Temp <- Temp * control.SA$coolingRate
               #Temp <- Temp * (1 + (F_c - F_b)/F_c )
             }
           } else {
             Temp <- control.SA$init.Temp
             while (Temp > control.SA$stop.Temp){
               M <- 1
               accepted.M <- 1
               while (M <= control.SA$Mtrial ){
                 which_var <- sample(c(1:nvar),1, prob=level.num)
                 first.all.neighbour <- partition_all_neighbour_restricted(x = S_c_list[[which_var]],
                                                                           group=group[[which_var]],
                                                                           apart=apart[[which_var]],
                                                                           method=transition.method,
                                                                           IncludeOrigin = IncludeOrigin)
                 first.neighbour.num <- sample(c(1:nrow(first.all.neighbour)), 1, replace=TRUE)
                 S_n_list <- S_c_list
                 S_n_list[[which_var]] <- first.all.neighbour[first.neighbour.num, ]
                 F_n <- fc.model.refit(varia.list=varia.list,
                                       merge.list=S_n_list,
                                       mod = model)[[1]]
                 All_states <- append(All_states, list(S_n_list))
                 All_bic <- c(All_bic, F_n)
                 if (F_n < F_c || stats::runif(1, 0, 1) < exp( (F_c - F_n) / Temp)) {
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
               #accepted.vec <- c(accepted.vec, accepted.M)
               Temp <- Temp * control.SA$coolingRate
               #Temp <- Temp * (1 + (F_c - F_b)/F_c )
             }
           }
           take.bic <- sort(unique(Accept_bic))
           take.weight <- BMAweight_bic( sort(unique(Accept_bic), decreasing=FALSE) )
           take.state <- NULL
           for (i in seq_len(length(take.bic))){
             temp.state <- unlist(Accept_states[ which( Accept_bic == sort(unique(Accept_bic))[i] )[1] ], recursive = FALSE)
             take.state <- append(take.state, list(temp.state))
           }
           Graycode.vec <- NULL
           Partition.vec <- NULL
           for (a in c(1:length(take.bic))){
             r1 <- NULL; r2 <- NULL
             mlist <- take.state[[a]]
             for (i in c(1:length(mlist))){
               r1 <- cbind(r1, mlist[i])
               r2 <- cbind(r2, graycode_to_partition(mlist[[i]]))
             }
             Graycode.vec <- rbind(Graycode.vec, r1)
             Partition.vec <- rbind(Partition.vec, r2)
           }
           Table <- as.data.frame(cbind(Graycode.vec,
                                        Partition.vec,
                                        as.numeric(take.bic),
                                        as.numeric(take.weight)))
           names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
           rownames(Table) <- NULL

           if (!AddOnSearch){
             result <- list(method = "SA",
                            variable = varia.list,
                            best.bic = F_b,
                            best.state = S_best_list,
                            all.bic = take.bic,
                            all.states = Graycode.vec,
                            all.weights = take.weight,
                            all.partitions = Partition.vec,
                            table = Table,
                            model = model,
                            addon = FALSE)
             #accepted.vec = accepted.vec
           } else{
             if (verbose) cat("Additional greedy search is used.", "\n")
             assign("addon.bic.trace", NULL, envir = environment(FCBMA))
             assign("addon.code.trace", NULL, envir = environment(FCBMA))
             res.add <- AddOnSearch(merge.code = S_best_list,
                                    varia.list = varia.list,
                                    model = model,
                                    transition.method=transition.method,
                                    IncludeOrigin = IncludeOrigin)
             addon.bestbic <- res.add[[1]]
             addon.beststate <- res.add[[2]]
             addon.bic.trace <- get("addon.bic.trace", envir = environment(FCBMA))
             addon.code.trace <- get("addon.code.trace", envir = environment(FCBMA))
             if (is.null(addon.bic.trace)){
               result <- list(method = "SA",
                              variable = varia.list,
                              best.bic = F_b,
                              best.state = S_best_list,
                              all.bic = take.bic,
                              all.states = Graycode.vec,
                              all.weights = take.weight,
                              all.partitions = Partition.vec,
                              table = Table,
                              model = model,
                              addon = TRUE)
             } else {
               #-------------------
               indx <- order(c(take.bic, addon.bic.trace), decreasing = FALSE)
               all.bic <- c(take.bic, addon.bic.trace)[indx]
               all.weights <- BMAweight_bic(all.bic)
               new.graycode.vec <- NULL
               new.partition.vec <- NULL
               for (a in c(1:length(addon.bic.trace))){
                 r1 <- NULL; r2 <- NULL
                 mlist <- addon.code.trace[[a]]
                 for (i in c(1:length(mlist))){
                   r1 <- cbind(r1, mlist[i])
                   r2 <- cbind(r2, graycode_to_partition(mlist[[i]]))
                 }
                 new.graycode.vec <- rbind(new.graycode.vec, r1)
                 new.partition.vec <- rbind(new.partition.vec, r2)
               }
               all.states <- rbind(Graycode.vec, new.graycode.vec)[indx,]
               all.partitions <- rbind(Partition.vec, new.partition.vec)[indx,]
               Table <- as.data.frame(cbind(all.states,
                                            all.partitions,
                                            as.numeric(all.bic),
                                            as.numeric(all.weights)))
               names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
               rownames(Table) <- NULL
               #------------------
               result <- list(method = "SA",
                              variable = varia.list,
                              best.bic = addon.bestbic,
                              best.state = addon.beststate,
                              all.bic = all.bic,
                              all.states = all.states,
                              all.weights = all.weights,
                              all.partitions = all.partitions,
                              table = Table,
                              model = model,
                              addon = TRUE)
             }
           }
           structure(result, class = c('FCBMA'))
         },
         GA = {
           newdat <- eval(stats::getCall(model)$data, envir = environment(stats::formula(model)))
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
             for (each in c(1:control.GA$popnSize)){
               S_c_list <- NULL
               for (j in seq_len(nvar)){
                 S_c_temp <- convert_canonical_graycode( sample(c(1:level.num[j]), size = level.num[j], replace=TRUE) )
                 S_c_list <- append(S_c_list, list(S_c_temp))
               }
               Popn <- append(Popn, list(S_c_list))
             }
           } else {
             for (each in c(1:control.GA$popnSize)){
               S_c_list <- NULL
               for (j in seq_len(nvar)){
                 svec <- sample(c(1:level.num[j]), size = level.num[j], replace=TRUE)
                 svec[group[[j]]] <- svec[group[[j]]][1]
                 svec[apart[[j]]] <- seq_len(length(apart[[j]]))
                 S_c_temp <- convert_canonical_graycode(svec)
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
             all.neighbour <- partition_all_neighbour(x = collapse.old,
                                                      method = transition.method,
                                                      IncludeOrigin = IncludeOrigin)
             collapse.new <- all.neighbour[ sample(c(1:nrow(all.neighbour)), 1) , ]
             thelist[[which.var]] <- collapse.new
             return(thelist)
           }
           Mutation.function.restricted <- function(thelist){
             which.var <- sample( c(1:nvar),1, prob=level.num )
             collapse.old <- unlist(thelist[[which.var]])
             all.neighbour <- partition_all_neighbour_restricted(x = collapse.old,
                                                                 group=group[[which.var]],
                                                                 apart=apart[[which.var]],
                                                                 method = transition.method,
                                                                 IncludeOrigin = IncludeOrigin)
             collapse.new <- all.neighbour[ sample(c(1:nrow(all.neighbour)), 1) , ]
             thelist[[which.var]] <- collapse.new
             return(thelist)
           }
           if (is.null(unique(unlist(group)))&&is.null(unique(unlist(apart)))){
             while (generation.k <= control.GA$MaxGen){
               # elitism
               elitism.num <- round(control.GA$popnSize*control.GA$elitism)
               NewPopn <- Popn[order(fitness.vec, decreasing=FALSE)[1:elitism.num]]
               Popn[ order(fitness.vec, decreasing=FALSE)[1:elitism.num] ] <-NULL
               OldPopn <- Popn
               OldFitness <- fitness.vec[-(order(fitness.vec, decreasing=FALSE)[1:elitism.num])]
               # copy # using tounament rule
               copy.popn <- OldPopn[TournamentSelection(round((1 - control.GA$CrossOverRate)*length(OldPopn)), k=2, OldFitness) ]
               # using roulette rule #copy.popn <- OldPopn[sample(c(1:length(OldFitness)), size = round((1- CrossOverRate)*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
               NewPopn <- append(NewPopn, copy.popn)
               # cross-over # one point cross over first using tournament rule
               cross.mat <- OldPopn[TournamentSelection(round(control.GA$CrossOverRate*length(OldPopn)), k=2, OldFitness)]
               # using roulette rule
               #cross.mat <- OldPopn[sample(c(1:length(OldFitness)), size = round(CrossOverRate*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
               crossover.popn <- Crossover.function(cross.mat)
               NewPopn <- append(NewPopn, crossover.popn)
               # mutation
               mutation.vec <- sample(c( (round(elitism.num)+1):(length(NewPopn)) ), round(control.GA$MutationRate*control.GA$popnSize), replace = FALSE)
               mutation.mat <- NewPopn[mutation.vec]
               NewPopn[mutation.vec] <- NULL
               mutation.popn <- lapply(mutation.mat, Mutation.function)
               NewPopn <- append(NewPopn, mutation.popn)
               fitness.vec <- sapply(NewPopn, fitness.calc.function)
               fittestOne <- c(fittestOne, min(fitness.vec))
               Popn <- NewPopn
               generation.k <- generation.k + 1
               if (verbose) {cat("generation", generation.k, "\n")}
             }
           } else {
             while (generation.k <= control.GA$MaxGen){
               # elitism
               elitism.num <- round(control.GA$popnSize*control.GA$elitism)
               NewPopn <- Popn[order(fitness.vec, decreasing=FALSE)[1:elitism.num]]
               Popn[ order(fitness.vec, decreasing=FALSE)[1:elitism.num] ] <- NULL
               OldPopn <- Popn
               OldFitness <- fitness.vec[-(order(fitness.vec, decreasing=FALSE)[1:elitism.num])]
               # copy # using tounament rule
               copy.popn <- OldPopn[TournamentSelection(round((1- control.GA$CrossOverRate)*length(OldPopn)), k=2, OldFitness) ]
               # using roulette rule #copy.popn <- OldPopn[sample(c(1:length(OldFitness)), size = round((1- CrossOverRate)*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
               NewPopn <- append(NewPopn, copy.popn)
               # cross-over # one point cross over first using tournament rule
               cross.mat <- OldPopn[TournamentSelection(round(control.GA$CrossOverRate*length(OldPopn)), k=2, OldFitness)]
               # using roulette rule
               #cross.mat <- OldPopn[sample(c(1:length(OldFitness)), size = round(CrossOverRate*dim(OldPopn)[1]) ,prob = OldFitness, replace = FALSE), ]
               crossover.popn <- Crossover.function(cross.mat)
               NewPopn <- append(NewPopn, crossover.popn)
               # mutation
               mutation.vec <- sample(c( (round(elitism.num)+1):(length(NewPopn)) ), round(control.GA$MutationRate*control.GA$popnSize), replace = FALSE)
               mutation.mat <- NewPopn[mutation.vec]
               NewPopn[mutation.vec] <- NULL
               mutation.popn <- lapply(mutation.mat, Mutation.function.restricted)
               NewPopn <- append(NewPopn, mutation.popn)
               fitness.vec <- sapply(NewPopn, fitness.calc.function)
               fittestOne <- c(fittestOne, min(fitness.vec))
               Popn <- NewPopn
               generation.k <- generation.k + 1
               if (verbose) {cat("generation", generation.k, "\n")}
             }
           }
           Best.State <- unlist(Popn[which.min(fitness.vec)], recursive = FALSE)
           Best.BIC <- min(fitness.vec)
           All.Popn.States <- Popn
           All.BIC <- fitness.vec
           Best.BIC.trace <- fittestOne
           take.bic <- sort(unique(All.BIC))
           take.weight <- BMAweight_bic( sort(unique(All.BIC), decreasing=FALSE) )
           take.state <- NULL
           for (i in seq_len(length(take.bic))){
             temp.state <- unlist(All.Popn.States[ which(All.BIC==sort(unique(All.BIC))[i] )[1] ], recursive = FALSE)
             take.state <- append(take.state, list(temp.state))
           }
           Graycode.vec <- NULL
           Partition.vec <- NULL
           for (a in c(1:length(take.bic))){
             r1 <- NULL; r2 <- NULL
             mlist <- take.state[[a]]
             for (i in c(1:length(mlist))){
               r1 <- cbind(r1, mlist[i])
               r2 <- cbind(r2, graycode_to_partition(mlist[[i]]))
             }
             Graycode.vec <- rbind(Graycode.vec, r1)
             Partition.vec <- rbind(Partition.vec, r2)
           }
           Table <- as.data.frame(cbind(Graycode.vec,
                                        Partition.vec,
                                        as.numeric(take.bic),
                                        as.numeric(take.weight)))
           names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
           rownames(Table) <- NULL

           if (!AddOnSearch){
             result <- list(method = "GA",
                            variable = varia.list,
                            best.bic = Best.BIC,
                            best.state = Best.State,
                            #all.bic = All.BIC,
                            #all.popn.states = All.Popn.States,
                            all.bic = take.bic,
                            all.states = Graycode.vec,
                            all.weights = take.weight,
                            all.partitions = Partition.vec,
                            table = Table,
                            model = model,
                            addon = FALSE)
           } else{
             if (verbose) cat("Additional greedy search is used.", "\n")
             assign("addon.bic.trace", NULL, envir = environment(FCBMA))
             assign("addon.code.trace", NULL, envir = environment(FCBMA))
             res.add <- AddOnSearch(merge.code = Best.State,
                                    varia.list = varia.list,
                                    model = model,
                                    transition.method=transition.method,
                                    IncludeOrigin = IncludeOrigin)
             addon.bestbic <- res.add[[1]]
             addon.beststate <- res.add[[2]]
             addon.bic.trace <- get("addon.bic.trace", envir = environment(FCBMA))
             addon.code.trace <- get("addon.code.trace", envir = environment(FCBMA))
             if (is.null(addon.bic.trace)){
               result <- list(method = "GA",
                              variable = varia.list,
                              best.bic = Best.BIC,
                              best.state = Best.State,
                              #all.bic = All.BIC,
                              #all.popn.states = All.Popn.States,
                              all.bic = take.bic,
                              all.states = Graycode.vec,
                              all.weights = take.weight,
                              all.partitions = Partition.vec,
                              table = Table,
                              model = model,
                              addon = FALSE)
             } else {
               #-------------------
               indx <- order(c(take.bic, addon.bic.trace), decreasing = FALSE)
               all.bic <- c(take.bic, addon.bic.trace)[indx]
               all.weights <- BMAweight_bic(all.bic)
               new.graycode.vec <- NULL
               new.partition.vec <- NULL
               for (a in c(1:length(addon.bic.trace))){
                 r1 <- NULL; r2 <- NULL
                 mlist <- addon.code.trace[[a]]
                 for (i in c(1:length(mlist))){
                   r1 <- cbind(r1, mlist[i])
                   r2 <- cbind(r2, graycode_to_partition(mlist[[i]]))
                 }
                 new.graycode.vec <- rbind(new.graycode.vec, r1)
                 new.partition.vec <- rbind(new.partition.vec, r2)
               }
               all.states <- rbind(Graycode.vec, new.graycode.vec)[indx,]
               all.partitions <- rbind(Partition.vec, new.partition.vec)[indx,]
               Table <- as.data.frame(cbind(all.states,
                                            all.partitions,
                                            as.numeric(all.bic),
                                            as.numeric(all.weights)))
               names(Table) <- c(paste(rep("Graycode",nvar),seq_len(nvar)), paste(rep("Partition", nvar),seq_len(nvar)), "BIC", "Model.weight")
               rownames(Table) <- NULL
               #------------------
               result <- list(method = "GA",
                              variable = varia.list,
                              best.bic = addon.bestbic,
                              best.state = addon.beststate,
                              all.bic = all.bic,
                              all.states = all.states,
                              all.weights = all.weights,
                              all.partitions = all.partitions,
                              table = Table,
                              model = model,
                              addon = TRUE)
             }
           }
           structure(result, class = c('FCBMA'))
         },
         stop("invalid method in FCBMA. Please use 'complete', 'SA' or 'GA'.") )
}



#' @export
print.FCBMA <- function(x, ...){
  txt <- paste0("'FCBMA' method with variable '",
                paste(x$variable, collapse="', '"), "' collapsed")
  cat(txt, "\n")
  cat("\n")
  cat("Available components:\n")
  print(names(x))
  invisible(x)
}
