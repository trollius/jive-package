#' @title Jive MCMC
#' @description Implements Markov chain Monte Carlo sampling for trait evolutionary models with intraspecific data
#' 
#' @details This function runs MCMC sampling on jive object \code{\link{jiveMake}}. The jive object contains 
#' both the dataset and set of model to be used in MCMC. This function implements both a conventional MCMC
#' and an MCMC with thermodynamic integration. The latter option is turned off by default and can be changed by
#' setting ncat to values > 1. The recommended ncat for TI is 10. When setting ncat > 1, make sure to specify burning.
#' As a rule of thumb set burning to 1/10 fraction of ngen. 
#' 
#' @param jive an object of class "jive" (see details)
#' @param log.file name of the output file that will store the log of MCMC chain
#' @param sampling.freq sampling frequency of the MCMC chain (how often chain will be saved into output file
#' @param print.freq printing frequency of the MCMC chain (how often chain will be printed in the R console)					
#' @param ncat number of classes for thermodynamic integration (see details)
#' @param beta.param beta value to define classes for thermodynamic integration (see details)
#' @param ngen number of generation in MCMC chain
#' @param burnin a burning phase of MCMC chain (has to be specified for thermodynamic integration)
#' @param update.freq update frequencies for likelihood and prior level parameters
#' @export
#' @author Anna Kostikova and Daniele Silvestro
#' @return none
#' @examples
#' ## Load test data
#' data(traitsOU1)
#' data(treeOU1)
#' ## Run a simple MCMC chain
#' my.jive <- jiveMake(treeOU1, traitsOU1,  model.var="OU1", model.mean="BM", model.lik="Multinorm")
#' jiveMCMC(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 
#'
#' ## Run an MCMC chain with thermodynamic integration
#' jiveMCMC(my.jive, log.file="my.jive_MCMC.log", ncat=10, sampling.freq=10, print.freq=100, ngen=5000, burnin=500) 


# MCMC MAIN ALGORITHM
jiveMCMC <- function(jive, log.file="jive_mcmc.log", sampling.freq=1000, print.freq=1000, 
				ncat=1, beta.param=0.3, ngen=5000000, burnin=0, update.freq=NULL)
{
	
			 
	# counter for the total chain length 
	real.iter <- 1 

	# here we define the chain length for each category
	it <- (ngen - burnin)/ncat 
		 if (ncat > 1) {
			  # here we get the heating parameter for a chain - scaling classes (MIGRATE)
			  beta.class <- defClasses(ncat, beta.param) 
		 } else {
			  beta.class <- 1
		 }
		 # start looping over scaling classes 
		 for (i.beta in 1:length(beta.class)){ 
			  # apply burning only to a first class
			  if (i.beta > 1) { 
				   burnin <- 0
			  }
			  #exctract temperature
			  temperature <- beta.class[i.beta]
				   acc <- 0 #acceptance count
				   for (iteration in 1:(it + burnin)) {
						hasting.ratio <- 0
						if (real.iter == 1){

							# initialize update frequencies
							update.freq <- initUpdateFreq(update.freq)

							
							# initialize MCMC parameters
							mspA  <- jive$lik$mspinit # initialize means for species
							sspA  <- jive$lik$sspinit # initialize sigma.sq for species
							smvnA <- jive$prior_mean$init[1] # could be a random number, initialize sigma.sq for MVN
							mmvnA <- jive$prior_mean$init[2] # could be a random number, initialize mean for MVN
							bmouA <- jive$prior_var$init # could be aither more realistic values such as means and sds of true data
							 
						}
						
						msp  <- mspA
						ssp  <- sspA
						mmvn <- mmvnA
						smvn <- smvnA
						bmou <- bmouA
						
						r    <- runif(1)
						
						# Update on all parameters level

						# update msp, ssp 
						if (r < update.freq[1]) {
						
							 # 5 is just for now number of species updated
							 ind <- sample(1:length(msp), 5, replace=F)
							
							 if (runif(1) < 0.5) {
								  # updating random 5 values from the vector of means 
								  msp[ind] <- slidingWin(mspA[ind], jive$lik$mspws[ind]) 
							 } else {
								  # updating random 5 values from the vector of sigmas
								  ssp[ind] <- abs(slidingWin(sspA[ind], jive$lik$sspws[ind]))
							 }
						}
						# update MVN parameters <-----------------------------  code refactoring is needed to allow other than BM models
						else if (r < update.freq[2]) {
						
							 if (runif(1) < 0.5){
								  # updating mmvn - negative values allowed
								  mmvn <- slidingWin(mmvnA, jive$prior_mean$ws[2]) 

							 } else {
								  # updating smvn
								  smvn <- abs(slidingWin(smvnA, jive$prior_mean$ws[1]) )
							 }
							 
						} else {# update BMOU parameters 
												 
							 ind <- sample(1:length(jive$prior_var$ws), 1)
							 bmou[ind] <- abs(slidingWin(bmouA[ind], jive$prior_var$ws[ind])) # updating bmou parameters
							 
						}
						

						# Hyperpriors on all parameters level
						# mean of MVN can be negative due to PCA - mean hprior

						hprior_mean <- mapply(do.call, jive$prior_mean$hprior, lapply(c(smvn, mmvn), list))
						hprior_var <- mapply(do.call, jive$prior_var$hprior, lapply(bmou, list))
						hprior <- c(hprior_mean, hprior_var)
						
						if (real.iter > 1) {
							 Lik <- LikA
							 Prior_mean <- Prior_meanA #priorMVNA
							 Prior_var <- Prior_varA #priorBMOUA
						}
						
						# do this for first step always (because we need to have all probabiliiles)     
						if (r < update.freq[1] || real.iter == 1) {
						
							 Lik 		<- jive$lik$model(msp, ssp, jive$data$traits, jive$data$counts) # traits, counts
							 #print(Lik)
							 Prior_mean <- jive$prior_mean$model(c(smvn, mmvn), msp, jive$data$tree) # tree Conditional prior level
							 #print(paste("Prior mean ", Prior_mean))
							 
							 Prior_var 	<- jive$prior_var$model(bmou, ssp, jive$data$tree, jive$data$map)
							 #print(paste("Prior var ", Prior_var))
							 										 
						} else if (r<update.freq[2]) {
							 Prior_mean <- jive$prior_mean$model(c(smvn, mmvn), msp, jive$data$tree)
						} else {
							 
							 Prior_var 	<- jive$prior_var$model(bmou, ssp, jive$data$tree, jive$data$map)
						}

						# Posterior calculation

						# jsut for 1 real.iter we need to copy all calculated likelihoods and priors
						if (real.iter == 1) {
							 LikA <- Lik
							 Prior_meanA <- Prior_mean
							 Prior_varA <- Prior_var
							 hpriorA <- hprior
							 postA <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))
						}
						post <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))
						# acceptance probability
						tryCatch(
						{
						if (post - postA + hasting.ratio >= log(runif(1))){
							 acc = acc + 1
							 LikA = Lik
							 Prior_meanA = Prior_mean
							 Prior_varA = Prior_var
							 hpriorA = hprior
							 postA = post
							 mspA = msp
							 sspA = ssp
							 mmvnA = mmvn
							 smvnA = smvn
							 bmouA = bmou
							}
						}
						,error = function(e) NULL
						)

						# log to file with frequency sampling.freq
						if (real.iter == 1){
							cat("generation",'\t',"posterior",'\n')
							cat(sprintf("%s\t", jive$prior_var$header), "\n", append=FALSE, file=log.file)
							} 
						
						if (real.iter %% sampling.freq == 0 & real.iter >= burnin) {
							 cat(sprintf("%s\t", c(real.iter, postA, sum(LikA),
							 Prior_meanA, Prior_varA, sum(hpriorA), smvnA, mmvnA, 
							 bmouA, mspA, sspA, (acc/iteration), temperature)),
							 "\n", append=TRUE, file=log.file) 
						}
						
						if (real.iter %% print.freq == 0 & real.iter >= burnin) {
							cat(real.iter,'\t',postA,'\n') 
						}
						  
				   real.iter = real.iter + 1
				   }
		 }# end of temperature
}



