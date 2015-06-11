function() {
	prof <- profileExample(n = 50, p = 250, n.burn = 100, n.reps = 200)
	rows <- grep("elementwise_DP.functions.R", row.names(prof))


	prof <- prof[order(expand_numbers(row.names(prof))),]


}


expand_numbers <- function(names) {
	names <- sub("R#(\\d)$", "R#000\\1", names)
	names <- sub("R#(\\d\\d)$", "R#00\\1", names)
	names <- sub("R#(\\d\\d\\d)$", "R#0\\1", names)
	names
}


quick_profile <- function() {

  set.seed(666)
  simulation <- simulateExample(n = 25, p = 125)
  system.time(
      posterior <- fitExample(simulation, n.burn = 100, n.reps = 200, computeMode = "R")
  )
  d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
  mean.taxicab <- mean(posterior$mean.taxicab.v)
  se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
  d_credible.v
  mean.taxicab
  se_mean.taxicab

  set.seed(666)
  simulation <- simulateExample(n = 25, p = 125)
  system.time(
    posterior <- fitExample(simulation, n.burn = 100, n.reps = 200, computeMode = "C")
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
