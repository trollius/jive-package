require(phytools)
require(MASS)
require(OUwie)
## number of species we want to simulate
n <- 50

## generate tree with a pure birth model and scale it to the height of 1
tree <- pbtree(b = 1, n = n, scale = 1, nsim = 1, ape = TRUE)

## set parameters for OU1 model of species-specific variances
sig.sq <- 0.9
alpha  <- 0.1
theta0 <- 1
theta  <- 5

## set parameters for BM model of specific-specific means
sig.sq.bm <- 0.5
mu0       <- 350

## set mean number of observations per species
mean.obs <- 20

## get selective regimes (all 1s because of OU1 model)
y <- data.frame(tree$tip.label, rep(1, n))

## add node labels
tree$node.label <- rep("1", n-1)

## simulate species-specific variances 
sigma.val <- abs(OUwie.sim(tree, y, simmap.tree=FALSE, scaleHeight=TRUE, alpha=rep(alpha,2), sigma.sq=rep(sig.sq,2), theta0=theta0, theta=theta)$X)

## simulate species-specific means
mean.val <- mvrnorm(mu=rep(mu0, length(tree$tip)), Sigma=(sig.sq.bm * vcv(tree)))

## draw a random number of intraspecific observations for each species
spec.obs <- rpois(n, mean.obs)

## generate a data matrix where rows are species and columns are individual observations	
traits <- matrix(rnorm(max(spec.obs) * n, mean=mean.val, sd=sqrt(sigma.val)), nrow=n, ncol=max(spec.obs))
traits <- cbind(as.matrix(max(spec.obs) - spec.obs), traits)

## function to replace empty cells with NA
foo <- function(x){
	to <- x[1]
	x[1:(to + 1)] <- NA
	return(x[-1])
}

## apply to data matrix	
traits <- as.matrix(t(apply(traits, 1, foo)))

## add species names to rownames
rownames(traits) <- tree$tip.label
my.jive <- jiveMake(tree, traits,  model.var="OU1", model.mean="BM", model.lik="Multinorm")