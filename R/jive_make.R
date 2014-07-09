#' @title Make jive object
#' @description Creates a jive object from intraspecific observations and species phylogeny
#' 
#' @details This function creates a jive object needed for \code{\link{jiveMCMC}} function.  
#' Trait values must be stored as a matrix, where lines are vectors of observations for each species, with NA for no data.
#' Rownames are species names. Phylogenetic tree must be provided as either simmap object (for models with multiple regimes)
#' or as a phylo object (for BM or OU1 models). Rownames and tip labels of a phylogenetic tree should match exactly. 
#' There are three models implemeted for estimation of species variances evolution - BM, OU1 and OUM. Evolution of 
#' species means is only implemented with BM model. Species-specific distribution are models as multivariate normal distribution
#' 
#' 
#' @param simmap an object of class "jive" (see details)
#' @param traits name of the output file that will store the log of MCMC chain
#' @param model.var sampling frequency of the MCMC chain (how often chain will be saved into output file
#' @param model.mean printing frequency of the MCMC chain (how often chain will be printed in the R console)					
#' @param model.lik number of classes for thermodynamic integration (see details)
#' @export
#' @examples
#' my.jive <- jiveMake(phy, traits,  model.var="OU1", model.mean="BM", model.lik="Multinorm") 


# library(geiger)
# jive structure
# jive$data			# $traits 
					# $counts 
					# $tree  
					# $vcv 
					# $map
					# $nreg					
				
# jive$lik			# $model
					# $wsmsp
					# $wsssp
					# $initmsp
					# $initssp
				
# jive$prior_mean 	# $model
					# $init
					# $ws
					# $hprior
						
# jive$prior_var	# $model
					# $init
					# $ws
					# $hprior

					
# traits - matrix, where lines are vector of observations for each species, with NA for no data.
# phy - simmap object (OUM) or phylo object (for BM or OU1)
jiveMake <- function(simmap, traits, model.var="OU1", model.mean="BM", model.lik="Multinorm"){

	jive <- list()
	
	if (name.check(simmap, traits) != "OK") {
	
		stop("Species do not match in tree and traits")
	
	} else {

		traits<-traits[simmap$tip.label,]
		
		if (model.var %in% c("BM", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime# td$phy$node.label <- rep("1", n - 1)
			jive$data$map <- matrix(rep(1, ((dim(traits)[1]) * 2) - 2))
		} else {
			jive$data$map <- relSim(simmap)$mapped.edge
		}
				
		jive$data$traits 					<- traits
		jive$data$counts 					<- apply(traits, 1, function (x) {sum( !is.na(x) )})
		jive$data$tree   					<- simmap
		jive$data$vcv    					<- vcv(simmap)
		jive$data$nreg   					<- dim(jive$data$map)[2]
		
		print(jive$data$nreg)
		
		
		if (model.lik == "Multinorm") {
			jive$lik$model 					<- likMultinorm
			jive$lik$mspws 					<- initWinSizeMVN(jive$data$traits)$msp
			jive$lik$sspws 					<- initWinSizeMVN(jive$data$traits)$ssp
			jive$lik$mspinit				<- initParamMVN(jive$data$traits)$mspA
			jive$lik$sspinit				<- initParamMVN(jive$data$traits)$sspA
				
		}
		
		if (model.mean == "BM" ) {
			jive$prior_mean$model 			<- likBM
			jive$prior_mean$init  			<- initParamMBM(jive$data$traits)  # check
			jive$prior_mean$ws	  			<- initWinSizeMBM(jive$data$traits)  # check
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5))		  # sigma
			jive$prior_mean$hprior$m		<- make.hpfun("Uniform", c(-10000,10000)) # anc.mean

		
		}
		
		if (model.var == "BM" ) {
			jive$prior_var$model 			<- likBM
			jive$prior_var$init  			<- initParamVBM(jive$data$traits) # check
			jive$prior_var$ws	  			<- initWinSizeVBM(jive$data$traits) # check
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_var$hprior$m			<- make.hpfun("Gamma", c(1.1,5)) # anc.mean
			jive$prior_var$header			<- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
												"mbm_sig.sq", "mbm_anc.st",  "vbm_sig.sq", "vbm_anc.st", 
												paste(rownames(jive$data$traits), "_m", sep=""),
												paste(rownames(jive$data$traits), "_v", sep=""),
												"acc", "temperature")
	
		}

		
		if (model.var == "OU1" || model.var == "OUM") {
			jive$prior_var$model 			<- likOU
			jive$prior_var$init  			<- initParamVOU(jive$data$traits, jive$data$nreg)  # check
			jive$prior_var$ws	  			<- initWinSizeVOU(jive$data$traits, jive$data$nreg)  # check
			jive$prior_var$hprior$a			<- make.hpfun("Gamma", c(1.1,5)) # alpha
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_var$hprior$m			<- make.hpfun("Gamma", c(1.1,5)) # anc.mean
			for (i in 1:jive$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_var$hprior[[ti]]	<- make.hpfun("Gamma", c(1.1,5))
			}
			jive$prior_var$header			<- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
												"mbm_sig.sq", "mbm_anc.st", "vou_alpha", "vou_sig.sq", "vou_anc.st", 
												paste("vou_theta", seq(1:jive$data$nreg),sep=""),
												paste(rownames(jive$data$traits), "_m", sep=""),
												paste(rownames(jive$data$traits), "_v", sep=""),
												"acc", "temperature")
		
		}
		
	}	
	return(jive)

}




