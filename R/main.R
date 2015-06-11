
#' @import MASS
#' @import mvtnorm
#' @export
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


#' Fit an example DPP model
#'
#' @description \code{fitExample} fits an example DPP model
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' simulation <- simulateExample(n = 25, p = 125)
#'
#' # Fit model
#' posterior <- fitExample(simulation, n.burn = 10, n.reps = 20)
#'
#' # Summarize posterior
#' d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
#' mean.taxicab <- mean(posterior$mean.taxicab.v)
#' se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
#' }
#'
#' @export
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

	posterior <- fn.mcmc(text="CLUST ANALYZE...",
											 data$X$true, data$X$data,
											 n.burn, n.reps, max.row.nbhd.size, row.frac.probes, col.frac.probes,
											 data$parm)
	return (posterior)
}

#' @export
profileExample <- function(n = 25,
													 p = 250,
													 n.burn = 10,
													 n.reps = 20) {
	simulation <- simulateExample(n, p)

	Rprof(line.profiling = TRUE, interval = 0.001)
	out <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps)
	Rprof(NULL)

	summaryRprof(lines = "show")$by.self
}

#' @useDynLib NPCluster, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
createEngine <- function(sort) {
	.createEngine(sort)
}

#' @export
accessEngine <- function(engine) {
  .accessEngine(engine)
}

#' @export
computePmfAndNeighborhoods <- function(engine, n0, vec.n, small) {
  .computePmfAndNeighborhoods(engine, n0, vec.n, small)
}
