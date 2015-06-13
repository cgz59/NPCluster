function() {
	prof <- profileExample(n = 400, p = 500, n.burn = 10, n.reps = 20,
	                       row.frac.probes = 0.25, col.frac.probes = 0.25,
	                       computeMode = "C")

	# rows <- grep("elementwise_DP.functions.R", row.names(prof))

	expand_numbers <- function(names) {
	  names <- sub("R#(\\d)$", "R#000\\1", names)
	  names <- sub("R#(\\d\\d)$", "R#00\\1", names)
	  names <- sub("R#(\\d\\d\\d)$", "R#0\\1", names)
	  names
	}

	filteredProf <- prof[prof$self.pct >= 1.00,]
	orderedProf <- filteredProf[order(expand_numbers(row.names(filteredProf))),]


}





quick_profile <- function() {

  n <- 100
  p <- 500
  n.burn <- 10
  n.reps <- 20

  n <- 25
  p <- 125

  set.seed(666)
  simulation <- simulateExample(n = n, p = p)
  system.time(
      posterior <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                              computeMode = createComputeMode())
  )
  d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
  mean.taxicab <- mean(posterior$mean.taxicab.v)
  se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
  d_credible.v
  mean.taxicab
  se_mean.taxicab

  set.seed(666)
  simulation <- simulateExample(n = n, p = p)
  system.time(
    posterior <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                            computeMode = createComputeMode(language = "C",
                                                            exactBitStream = TRUE))

  )
  d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
  mean.taxicab <- mean(posterior$mean.taxicab.v)
  se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
  d_credible.v
  mean.taxicab
  se_mean.taxicab

  set.seed(666)
  simulation <- simulateExample(n = n, p = p)
  system.time(
    posterior <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                            computeMode = createComputeMode(language = "C",
                                                            exactBitStream = FALSE))

  )
  d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
  mean.taxicab <- mean(posterior$mean.taxicab.v)
  se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
  d_credible.v
  mean.taxicab
  se_mean.taxicab





#   user  system elapsed
#   23.454   0.354  22.890
#   > d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
#   > mean.taxicab <- mean(posterior$mean.taxicab.v)
#   > se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
#   > d_credible.v
#   2.5%    97.5%
#     0.000000 0.318525
#   > mean.taxicab
#   [1] 0.0057728
#   > se_mean.taxicab
#   [1] 1.814736e-06

}
