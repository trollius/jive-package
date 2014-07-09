#' @title Simulate trees
#' @description What it does
#' 
#' @details This function 
#' 
#' 
#' 
#' 
#' @param nb_diversification_shift fff
#' @param nb_mass_extinction fff
#' @param total_nb_species fff
#' @param lambda_max ff		
#' @param lambda_min fff
#' @param mu_min ff
#' @param survival_rate_min ff
#' @param seed ff
#' @export
#' @examples
#' lol <- classical_model(1, 1, 500, 0.1, 0.01, .01, .2, seed=1506621)


classicalModel <- function(nb_diversification_shift, nb_mass_extinction, total_nb_species, lambda_max, lambda_min, mu_min, survival_rate_min, seed=1){
  set.seed(seed)
  print(seed)
  vect_nb_spec <- rep(0, nb_diversification_shift+1)
  while (sum(vect_nb_spec<2) | sum(vect_nb_spec)!=total_nb_species)   vect_nb_spec <- round(rdirichlet(1, rep(1, nb_diversification_shift+1))*total_nb_species)
  param <- data.frame(l=rep(0, nb_diversification_shift+1), m=rep(0, nb_diversification_shift+1), spec=t(vect_nb_spec))
  for (i in 1:(nb_diversification_shift+1)){
    param[i, "l"] <- runif(1, min=lambda_min, max=lambda_max)
    param[i, "m"] <- runif(1, min=mu_min, max=param[i, "l"])
  }
#  param <- data.frame(l=runif(nb_diversification_shift+1, min=lambda_min, max=lambda_max))
 # print(param)
 # param <- cbind(param, data.frame(m=runif(nb_diversification_shift+1, min=mu_min, max=param[1, "l"]), spec=t(vect_nb_spec)))
  print(param)
             
  rownames(param) <- paste("s", 1:(nb_diversification_shift+1), sep="")
  mee_times <- c(0, runif(nb_mass_extinction, min=0, max=min(log(param[,"spec"])/(param[,"l"]-param[,"m"]))))  
  surv_rate <- c(1, runif(nb_mass_extinction, min=survival_rate_min, max=1))
  mee_times <- sort(mee_times)
  spec_names <- c(0, cumsum(vect_nb_spec)) ##definition of parameters
  br.ti <- c(0)
  print(param)
  print(mee_times)
  print(surv_rate)
  while (max(br.ti) < tail(mee_times, 1) | !max(br.ti)){
    tree <- list()
    for (i in 1:(nb_diversification_shift+1)){ #simulation of all sub trees
      tmp <- sim.rateshift.taxa(param[i, "spec"], 1, rep(param[i, "l",], length(surv_rate)), rep(param[i, "m"], length(surv_rate)), surv_rate, mee_times, complete = TRUE, K=0, norm = TRUE)
      tree[[i]] <- list()
      if (Ntip(tmp[[1]][[1]]) > vect_nb_spec[i]) tree[[i]][[1]] <- drop.tip(tmp[[1]][[1]], (vect_nb_spec[i]+1):Ntip(tmp[[1]][[1]])) #deletes extinct tips
      else tree[[i]][[1]] <- tmp[[1]][[1]]
      tree[[i]][[1]]$tip.label <- paste("t",(spec_names[i]+1):spec_names[i+1], sep="")
      tree[[i]][[1]]$node.label <-  paste("n", (spec_names[i]+1):(spec_names[i+1]-1), sep="")
      tree[[i]][[1]]$node.label[1] <- paste("s", i, sep="") #we mark the root as special to keep track of the changes in lambda mu
      tree[[i]][[2]] <- tmp[[2]] #we keep the root length to bind the trees in the next step
    }
    br.ti <- sapply(tree, function(x) branching.times(x[[1]])[1])
  }
  root_times <- sapply(tree, function(x) x[[2]])
  param <- cbind(param, root_times, mrca_times=br.ti)
  print(param)
  return(c(sorting_and_binding_all_trees(tree, param), list(mee=mee_times, seed=seed, survival=surv_rate)))
}



