
#' @import MASS
#' @import mvtnorm
#' @export
simulateExample <- function(n = 25, p = 250, prop.X.miss=0, tau = 0.5, tau_0 = 1.25) {

	###################
	# generate covariates adding random noise of specified level
	# create objects data and true
	###################

	true_parm <- gen.clust(n, p)

	true_parm$tau <- tau
	true_parm$tau_0 <- tau_0

	sim.X <- gen.X(n, p, prop.X.miss, true_parm)

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
#' @useDynLib NPCluster, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @export
fitExample <- function(data,
											 n.burn = 10,
											 n.reps = 20,
											 max.row.nbhd.size = 50, # should be small compared to n2 * G
											 max.col.nbhd.size = 25, # should be small compared to p
											 row.frac.probes = 0.05,
											 col.frac.probes = .1,
                       prob.compute.col.nbhd=.2,
											 dahl.flag=FALSE,
											 standardize.X=FALSE,
											 tBB_flag=FALSE,
											 computeMode = createComputeMode()) {

	if (!inherits(data, "NPClustSimulation")) {
		stop("Wrong data structure")
	}

  if (!inherits(computeMode, "computeMode")) {
    stop("Wrong compute mode")
  }

	###################
	# Detect clusters
	###################

	posterior <- fn.mcmc(text="CLUST ANALYZE...",
											 data$X$true, data$X$data,
											 n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes,
											 prob.compute.col.nbhd, data$parm, dahl.flag=dahl.flag, standardize.X, tBB_flag, computeMode)
	return (posterior)
}

#' @export
profileExample <- function(n = 25,
													 p = 250,
													 n.burn = 10,
													 n.reps = 20,
													 row.frac.probes = 0.05,
													 col.frac.probes = 0.05,
													 computeMode = createComputeMode(),
													 filename = "Rprof.out") {
	simulation <- simulateExample(n, p)

	Rprof(filename = filename, line.profiling = TRUE, interval = 0.001)
	posterior <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
	           row.frac.probes = row.frac.probes,
	           col.frac.probes = col.frac.probes,
	           computeMode = computeMode)
	Rprof(NULL)
	#summaryRprof(lines = "show")$by.self
	return(posterior)
}


#' createComputeMode
#' @export
createComputeMode <- function(language = "R",
                              exactBitStream = FALSE,
                              extraSort = TRUE,
                              completeTest = FALSE,
                              tolerance = 1E-10,
                              test1 = FALSE,
                              test2 = FALSE,
                              test3 = FALSE) {
  if (!(language %in% c("C","R"))) {
    stop("Invalid language")
  }

  useR <- (language == "R")
  device <- NULL
  if (!useR) {
    doSort <- (exactBitStream | extraSort)
    device <- .createEngine(doSort)
  }

  object <- list(
    computeR = (language == "R" | completeTest),
    computeC = (language == "C"),
    device = device,
    exactBitStream = exactBitStream,
    extraSort = extraSort,
    tolerance = tolerance,
    test1 = test1,
    test2 = test2,
    test3 = test3
  )
  class(object) <- "computeMode"
  return(object)
}

#' assertEqual
assertEqual <- function(x, y, tolerance = 0) {
  if (length(x) != length(y)) {
    stop(cat("C++ error -- length:", length(x), length(y)))
  }
  if (any(abs(x - y) > tolerance)) {
    stop(cat("C++ error -- value:", x, y, tolerance, sep = "\n"))
  }
}


