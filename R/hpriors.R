#' @title Title
#' @description Short description
#' 
#' @details Detailed description
#' 
#' 
#' 
#' 
#' @param hpf description
#' @param hp.pars description
#' @param ... description
#' @examples
#' my.hp <- make.hpfun(hpf="Uniform", hp.pars=c(1,2))

make.hpfun <-function(hpf="Uniform", hp.pars, ...){
	# Function that makes function of hyper prior
	#
	# Args:
	# 	hpf:  		distribution of hyper prior
	#	hp.pars:	parameters of a distribution
	#
	# Returns:
	#	Hyper-prior function (function).

		if (hpf == "Uniform"){
			my.f <- function(x, ...){
				hp <- sum(dunif(x, min=hp.pars[1], max=hp.pars[2], log=TRUE))
				return(hp)
			}
			
		}
		
		if (hpf == "Gamma"){
			my.f <- function(x, ...){
				hp <- sum(dgamma(x, shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
				return(hp)
			}
		}
		
		if (hpf == "Normal"){
			my.f <- function(x, ...){
				hp <- sum(dnorm(x, mean=hp.pars[1], sd=hp.pars[2], log=TRUE))
				return(hp)
			}
		}	
		
		
	return(my.f)
}
