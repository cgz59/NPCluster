
#' @import MASS
#' @import mvtnorm
simulateExample <- function(n = 25, p = 250, tau = 0.5, tau_0 = 1.25) {
		
	###################
	# generate covariates adding random noise of specified level
	# create objects data and true
	###################
	
	true_parm <- gen.clust(n, p)
	
	true_parm$tau <- tau
	true_parm$tau_0 <- tau_0
	
	sim.X <- gen.X(n, p, true_parm)

	simulation <- list(X = sim.X, parm = true_parm)
	class(simulation) <- "NPClustSimulation"
	
	return(simulation)
}


fitExample <- function(data, 
											 n.burn = 10, 
											 n.reps = 20,
											 max.row.nbhd.size = 50, # should be small compared to n2 * G
											 row.frac.probes = 0.05,
											 col.frac.probes = 0.05) {
	
	if (!inherits(data, "NPClustSimulation")) {
		stop("Wrong data structure")
	}

	###################
	# Detect clusters
	###################
	
	All.Stuff <- fn.mcmc(text="CLUST ANALYZE...",							
											 data$X$true, data$X$data,
											 n.burn, n.reps, max.row.nbhd.size, row.frac.probes, col.frac.probes, 
											 data$parm)
	
# 	d_credible.v <- quantile(All.Stuff$d.v, prob=c(.025,.975))
# 	
# 	mean.taxicab <- mean(All.Stuff$mean.taxicab.v)
# 	se_mean.taxicab <- sd(All.Stuff$mean.taxicab.v)/sqrt(n.reps)

	return (All.Stuff)
}


test <- function() {
	
	simulation <- simulateExample()
	posterior <- fitExample(simulation)
	
}