simple_check <- function(n = 25, p = 125,
                         n.burn = 100, n.reps = 500,
                         seed = 666, ...) {

  # n <- 25
  # p <- 125
  # n.burn <- 100
  # n.reps <- 500
  # seed <- 666

  # Run R only code
  set.seed(seed)
  simulation <- simulateExample(n = n, p = p, ...)
  t1 <- system.time(
      posterior1 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                               computeMode = createComputeMode())
  )
  G1 <- mean(posterior1$G.v)
  r1 <- mean(posterior1$rng)
  print(c(G1, r1))

  # Run R + C, but use some R values due to non-transitivity of floating-point ops
  set.seed(seed)
  simulation <- simulateExample(n = n, p = p, ...)
  mode <- createComputeMode(language = "C",
                            completeTest = TRUE,
                            exactBitStream = TRUE, # sets useCPdp <- FALSE
                            tolerance = 1E-10)
  t2 <- system.time(
    posterior2 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                             computeMode = mode)
  )
  G2 <- mean(posterior2$G.v) # Should exactly match mean(posterior1$G.v) # 53.886
  r2 <- mean(posterior2$rng) # Should exactly match mean(posterior1$rng)
  print(c(G2, r2))

  # Run R + C and check differences
  set.seed(seed)
  simulation <- simulateExample(n = n, p = p, ...)
  mode <- createComputeMode(language = "C",
                            completeTest = TRUE,
                            tolerance = 1E-10)

  t3 <- system.time(
    posterior3 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                             computeMode = mode)
  )
  G3 <- mean(posterior3$G.v)
  r3 <- mean(posterior3$rng)
  print(c(G3, r3))

  # Run C only
  set.seed(seed)
  simulation <- simulateExample(n = n, p = p, ...)
  mode <- createComputeMode(language = "C")

  t4 <- system.time(
    posterior4 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                             computeMode = mode)
  )
  G4 <- mean(posterior4$G.v)
  r4 <- mean(posterior4$rng)
  print(c(G4, r4))

  # Run multicore C
  set.seed(seed)
  simulation <- simulateExample(n = n, p = p, ...)
  mode <- createComputeMode(language = "C",
                            specialMode = "tbb")

  t5 <- system.time(
    posterior5 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
                             computeMode = mode)
  )
  G5 <- mean(posterior5$G.v) # Should exactly match mean(posterior4$G.v)
  r5 <- mean(posterior5$rng) # Should exactly match mean(posterior4$rng)
  print(c(G5, r5))

  # Run vectorized C   TODO -- THIS IS CURRENTLY BROKEN
  # set.seed(seed)
  # simulation <- simulateExample(n = n, p = p, ...)
  # mode <- createComputeMode(language = "C",
  #                           specialMode = "sse")
  #
  # t6 <- system.time(
  #   posterior6 <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
  #                            computeMode = mode)
  # )
  # G6 <- mean(posterior6$G.v) # Should exactly match mean(posterior4$G.v)
  # r6 <- mean(posterior6$rng) # Should exactly match mean(posterior4$rng)
  # print(c(G6, r6))

  t <- list(t1, t2, t3, t4, t5)
  G <- list(G1, G2, G3, G4, G5)
  r <- list(r1, r2, r3, r4, r5)

  list(t = t, G = G, r = r)
}

# Example execution
#
# NPCluster:::simple_check()
# NPCluster:::simple_check(n = 50)
# NPCluster:::simple_check(p = 250)
