
#' @import MASS
#' @import stats
#' @import mvtnorm
test <- function() {
	
# 	rm(list=ls())
	
# 	library(MASS)
# 	library(stats) 
# 	library(mvtnorm)
	
	options(error=recover)
	options(warn=2)
	
	# check before running
	n <- 50/2
	p <- 250/2
	
	###################
	# generate covariates adding random noise of specified level
	# create objects data and true
	###################
	
	true_parm <- gen.clust(n, p)
	#source("gen.clust.R")
	
	true_parm$tau <- .5
	true_parm$tau_0 <- 1.25
	
	sim.X <- gen.X(n, p, true_parm)
# 	source("gen.X.R")
	
	###################
	# Detect clusters
	###################
	
	n.burn <- 2500/250
	n.reps <- 5000/250
	
# 	source("iterations.r")
# 	source("elementwise_main.R")
# 	source("elementwise_DP.functions.R")
# 	source("PDP.functions.R")
 	
	max.row.nbhd.size <- 50 # should be small compared to n2*G
	row.frac.probes <- .05
	col.frac.probes <- .05
	
	All.Stuff <- fn.mcmc(text="CLUST ANALYZE...",
											 sim.X$true, sim.X$data, 
											 n.burn, n.reps, max.row.nbhd.size, row.frac.probes, col.frac.probes, 
											 true_parm)
	
	d_credible.v <- quantile(All.Stuff$d.v, prob=c(.025,.975))
	
	mean.taxicab <- mean(All.Stuff$mean.taxicab.v)
	se_mean.taxicab <- sd(All.Stuff$mean.taxicab.v)/sqrt(n.reps)

	return (list(
		parm = true_parm,
		data = sim.X$data,
		posterior = All.Stuff
		))
}
