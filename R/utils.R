

defClasses <- function(ncat=10, beta.param=0.3){ 
	# Defines classes for thermodynamic integration.
	# For details of the method see Xia et al 2011 Sys Bio.
	#
	# Args:
	# 	ncategories: number of classes that will be used in thermodynamic integration.
	#	beta.param:  parameter describing the shape of a beta distribution.
	# 
	# Returns:
	#	The vector of temperatures for thermodynamic integration.
	
	
    K <- ncat-1
    k <- 0:K
    b <- k/K
    temp<- rev(b^(1/beta.param))
	return(temp)
}

       
initUpdateFreq <- function(update.freq=NULL){
	# Initializes update frequencies for likelihood and two prior levels.
	#
	# Args:
	# 	update.freq: the vector (length = 3) of update frequencies (likelihood, priorMBM, priorVOU/VBM).
	#
	# Returns:
	#	The vector of update frequencies which sums to 1. 
	
	if (length(update.freq) != 3 && !is.null(update.freq)) {
		stop("Update.freq must contain 3 elements" )
	}
	
	
	# calculate update frequencies
	if (!is.null(update.freq)) {
		update.freq	<- cumsum(update.freq/sum(update.freq))
	} else {
		update.freq	<- cumsum(c(0.4,0.1,0.5))
	}
	
	return(update.freq)

}


initWinSize <- function(jive, window.sizes=NULL){
	# Calculates window sizes for sampling proposals.
	# User-defined window sizes are not supported at this stage.
	#
	# Args:
	# 	jive: jive.object (see makeJive function)
	#
	# Returns:
	#   A list of windows sizes for: msp - species-specific means, 
	#							     ssp - species-specific variances,
	#								 mvn - prior on means (MBM),
	#								 svn - prior on variances (VBM/VOU).
	
	ll.ws <- list()
	xx    <- apply(jive$traits, 1, sd, na.rm = T) # CHECK IF ITS SD OR VAR
	yy    <- sd(xx)
	
	ll.ws[["msp"]] <- xx 
	ll.ws[["ssp"]] <- xx
	ll.ws[["mvn"]] <- c(2 * yy, yy) # anc.state of means windows size, evol rate of means window size
	
	if (jive$model == "BM") {
		ll.ws[["svn"]] <- c(2 * yy, yy) # anc.state of sigmas windows size, evol rate of sgimas window size
	} else { 
		ll.ws[["svn"]] <- c(1.5, yy, 2 * yy, rep(2 * yy, jive$nreg)) # alpha from max likelihood on observed std dev <------------------------ alpha parameter to adjust
	}
	
	return(ll.ws)
	
}


## tranform simmap into map, input - simmap object
relSim <- function(x) {
	
	foo<-function(x) {
		x/sum(x)
	}
	
	x$mapped.edge <- t(apply(x$mapped.edge, 1, FUN=foo))
	x$mapped.edge <- x$mapped.edge[, order(colnames(x$mapped.edge))]
	
	return(x)
				
}


##-------------------------- initialize windows sizes functions
initWinSizeMVN <- function (x){

	ws		<- list()
	ws$msp	<- (apply(x, 1, sd, na.rm = T))/4 # <--------------- may need further tuning
	ws$ssp	<- apply(x, 1, sd, na.rm = T) 
	
	return(ws)

}
 
# input is trait matrix, rows are species, cols - observations
initWinSizeMBM <- function(x){
	
	xx <- sd(apply(x, 1, sd, na.rm = T)) # CHECK IF ITS SD OR VAR
	ws <- c(xx, 2 * xx) # evol rate of sigmas window size, anc.state of sigmas windows size, 

	return(ws)
	
}

initWinSizeVBM <- function(x){
	
	xx <- sd(apply(x, 1, sd, na.rm = T)) # CHECK IF ITS SD OR VAR
	ws <- c(xx, 2 * xx) # evol rate of sigmas window size, anc.state of sigmas windows size,

	return(ws)
	
}

# input is trait matrix, rows are species, cols - observations
initWinSizeVOU <- function(x, nreg){
	
	xx <- sd(apply(x, 1, sd, na.rm = T)) # CHECK IF ITS SD OR VAR
	ws <- c(1.5, xx, 2 * xx, rep(2 * xx, nreg)) # alpha from max likelihood on observed std dev <------------------------ alpha parameter to adjust
	
	return(ws)
	
}

##-------------------------- initialize start parameters values functions
# initialize MCMC parameters
initParamMVN <- function (x){

	init  <- list()
	init$mspA  <- apply(x, 1, mean, na.rm = T) # initialize means for species
	init$sspA  <- apply(x, 1, var, na.rm = T) # initialize sigma.sq for species
	
	return(init)
}


# initialize MCMC parameters				   
initParamMBM <- function(x){
		
	# initialize MCMC parameters
	init <- c(var(apply(x, 1, mean, na.rm = T)), mean(apply(x, 1, mean, na.rm = T))) 
	
	return(init)

}
				   
initParamVBM <- function(x){
		
	init <- c(runif(2, 0.5, 3)) # could be aither more realistic values such as means and sds of true data (
	#init <- c(2.941516,2.139533,1.299683,1.364224) just a check
	
	return(init)

}	

# initialize MCMC parameters (order: alpha, anc.state, sig, theta1, theta2...
initParamVOU <- function(x, nreg){
		
	init <- c(runif((nreg+3), 0.5, 3)) # could be aither more realistic values such as means and sds of true data
	#init <- c(2.941516,2.139533,1.299683,1.364224) just a check
	return(init)

}	



