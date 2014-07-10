#' @title Plot jive MCMC object
#' @description Plots jive MCMC output
#' 
#' @details This function plots estimated species-specific variances on the tree. The spec 
#' The species-specific variances are calculated by \code{\link{jiveProc}}
#' 
#' 
#' @param tree an object of class "jive" (see details)
#' @param proc.jive an object from \code{\link{jiveProc}} function (see details)
#' @param regime a named vector of selective regimes 
#' @param cols a vector of colours for selective regimes 					
#' @param cex.label magnification for tip labels
#' @param cex.circle magnification for circles representing species-specific variances
#' @param lab.off offset for tip labels
#' @param ladder if tree should be ladderized
#' @param ... additional parameters passed to \code{\link{ladderize}} and \code{\link{plot.phylo}} functions
#' @export
#' @author Anna Kostikova
#' @return nothing
#' @examples
#' jivePlot(phy1, my.l, regime) 



jivePlot <- function(tree, proc.jive, regime, cols=c("blue", "green"), cex.label=0.7, cex.circle=2, lab.off = 0.025, ladder=TRUE, ...){

		if (ladder){
			tree = ladderize(tree, right = TRUE)
		}
		regime <- regime[match(tree$tip.label, names(regime))]
		traits  <- proc.jive$var_pars[match(tree$tip.label, names(proc.jive$var_pars))]
		
		plot(tree,  label.offset = lab.off, show.tip.label = TRUE, cex = cex.label)
		tiplabels(pch = 19, col = cols[regime], cex = cex.circle * (traits/max(traits)))

}

