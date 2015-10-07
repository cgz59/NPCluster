

fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
	{

	# Computes multivariate normal density function
	# a little faster than dmvnorm function of R!

	if (missing(inv.sigma))
		{inv.sigma <- solve(sigma)
		}

	logdet <- as.numeric(determinant(inv.sigma, logarithm=TRUE)$mod)
	r <- length(x)
	Q <- colSums(inv.sigma * (x-mean))
	Q <- sum(Q * (x-mean))

	val <- -r/2*log(2*pi) + logdet/2 - Q/2
	if (!log)
		{val <- exp(val)
		}

	val
	}



fn.quality.check <- function(parm)
	{err <- 0

	if (!is.null(parm$clust$col.nbhd))
		{#if (sum(unlist(lapply(parm$clust$col.nbhd, length))) != length(parm$col.subset.I))
			#{err <- 1
			#}
		}

	if (sum(diag(parm$clust$tBB.mt) < 0) > 0)
		{err <- 2
		}

	if (ncol(parm$clust$A.mt) != parm$clust$G)
		{err <- 3
		}

	if (ncol(parm$clust$B.mt) != (parm$clust$G+1))
		{err <- 4
		}

	if ((sum(parm$clust$C.m.vec) + parm$clust$C.m0) != parm$p)
		{err <- 5
		}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 6
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 7
		}

	if ((sum(parm$clust$n.vec) + parm$clust$n0) != parm$N)
		{err <- 8
		}

	err
	}

######################

fn.init.clusters <- function(parm)
{

	X.mt <- parm$X
	num.centers <- parm$G.new

	options(warn=0)
	tmp2 <- kmeans(t(X.mt), iter.max=1000, centers=num.centers, nstart=2)
	options(warn=2)
	#
	parm$clust$c.v <- tmp2$cluster
	parm$clust$G <- length(tmp2$size)


	parm$clust$C.m.vec <- array(,parm$clust$G)

	for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 parm$clust$C.m.vec[g] <- sum(I.g)
		}

	###########################

	# start from PDP model with d=0 (i.e DP)
	parm$d <- 0


	parm
}


fn.eda <- function(parm, data, computeMode)
{

	parm <- fn.init.clusters(parm)
	# reintroduced on 6/29/12
	parm$G.max <- min(parm$p/2, round(parm$clust$G*1.1))
	parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)

	parm$Y <- parm$clust$A.mt <- array(,c(parm$n2,parm$clust$G))
	parm$clust$C.m.vec <- array(,parm$clust$G)

	for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 parm$clust$C.m.vec[g] <- m.g <- sum(I.g)
		x.g.v <- parm$X[,I.g]
		 if (m.g > 1)
			{x.g.v <- rowMeans(x.g.v)
			}
		parm$Y[,g] <- x.g.v
		}

	parm$clust$C.m0 <- parm$p - sum(parm$clust$C.m.vec)

	parm$clust$M <- parm$a.R
	parm$clust$M0 <- .01*parm$clust$M

	parm$clust$K <- data$K.max

	########## Fix  DP hyperprior

	parm$clust$mu2 <- mean(as.vector(parm$X))
	parm$clust$tau2 <- diff(range(as.vector(parm$X)))/6

	#################################

	parm$g <- rep(1:parm$clust$G,each=parm$n2)
	parm$N <- parm$clust$G * parm$n2

	parm$Y <- as.vector(parm$Y)

	# if (computeMode$useR) {

	tmp <- tryCatch({
	  iter.max <- ifelse((length(parm$Y) > 1000), 10, 1000)
	  kmeans(parm$Y, iter.max = iter.max, centers = data$K.max, nstart = 10,
	         algorithm = "Hartigan-Wong" # TODO: MAcQueen works better?
	  )}, error = function(e) {
	    print("Kmeans did not converge ... using random assignment")
	    cluster <- sample(1:data$K.max, size = length(parm$Y), replace = TRUE)
	    centers <- sapply(1:data$K.max, FUN = function(x) {
	      mean(parm$Y[which(cluster == x)])
	    })
	    list(cluster = cluster, centers = centers, size = data$K.max)
	  })

# 	} else {
#
# 	  tmp <- .fastKMeans(matrix(parm$Y, nrow = 1), data$K.max)
# 	  tmp$cluster <- tmp$cluster + 1
#
# 	}

	parm$clust$s.v <- tmp$cluster
	parm$clust$phi.v <- as.vector(tmp$centers)

	parm$clust$n.vec <- tmp$size
	# number of s equal to 0
	parm$clust$n0 <- 0

	parm$clust$s.mt <- array(parm$clust$s.v, c(parm$n2,parm$clust$G))

	for (g in 1:parm$clust$G)
		{parm$clust$A.mt[,g] <- parm$clust$phi.v[parm$clust$s.mt[,g]]
		}

	sum.resid.sq <- 0

	for (g in 1:parm$clust$G)
		{flag.v <- parm$clust$c.v == g
		X.g.mt <- parm$X[,flag.v]
		a.g.v <- parm$clust$A.mt[,g]
		resid.g.mt <- X.g.mt - a.g.v
		sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
		}

	parm$tau_int <- parm$tau <- sqrt(sum.resid.sq/parm$n2/parm$p)

	###################################

	parm$tau_0 <- sqrt(1+parm$tau^2)

	# 1-parm$tau^2/var(as.vector(parm$X))

	## objects of full size (based on all n2 cases)
	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
	parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt

	parm <- fn.assign.priors(parm, data)

	list(NA, parm)

	}



fn.gen.clust <- function(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)
	{

	# first, impute missing X values by their column-specific means
	# + a small error term to guarantee non-tied values
 	parm$X <- data$X
 	tmp.mean.v <- apply(data$X, 2, median, na.rm=TRUE)
	tmp.sd.v <- apply(data$X, 2, sd, na.rm=TRUE)
 	for (j in 1:parm$p)
		{indx.j <- is.na(data$X[,j])
		if (sum(indx.j) > 0)
			{parm$X[indx.j,j] <- tmp.mean.v[j] + rnorm(n=sum(indx.j), sd=tmp.sd.v[j]/5)
			}
		}

	parm$G.new <- data$G.max
	tmp <- fn.eda(parm, data, computeMode)
	parm <- tmp[[2]]

	#################

	parm <- fn.hyperparameters(data, parm)

	parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes=1, computeMode)

	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
	parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt

	parm

	}

fn.init <- function(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, computeMode = "R")
	{

#	parm <- true_parm

	parm <- NULL

	parm$n2 <- dim(data$X)[1] # TODO Check
	parm$p <- dim(data$X)[2]  # TODO Check

	# mass parameter of elementwise(s) groups
	# stored later in parm$clust$M
	parm$a.R <- true$a.R

	# mass paramater of columns
	parm$b1 <- true$b1

	# mass paramater of column-intercept cluster
	parm$b0 <- true$b0

#	parm$clust$C.m0 <- 0
#
#	parm$clust$A.mt <- array(,c(parm$n2,parm$clust$G))
#
#	for (g in 1:parm$clust$G)
#		{z.v <- parm$clust$s.mt[,g]>0
#		parm$clust$A.mt[z.v,g] <- parm$clust$phi.v[parm$clust$s.mt[z.v,g]]
#		parm$clust$A.mt[-z.v,g] <- 0
#		}
#
#	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
#
#	parm$shift <- true$shift
#
#	parm$tau_int <- parm$tau
#
#	parm$X <- data$X

	############################
	# For delta neighborhoods
	############################

	parm$col.delta <- .1

	# delta-neighborhood threshold for elements
	parm$row.delta <- .1

	#########################################
	# generating the R- and C- clusters
	########################################

	parm$shift <- true$shift
	parm <- fn.gen.clust(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)

	parm <- fn.assign.priors(parm, data)

	parm

	}


fn.assign.priors <- function(parm, data)
	{

	parm$prior$tau <- NULL
	parm$prior$tau$alpha.tau <- 1e-2
	parm$prior$tau$beta.tau <- 1e-2

	parm$prior$tau$max <- sqrt(.75)*sd(as.vector(data$X), na.rm=TRUE)
	parm$prior$tau$min <- 1e-10
	parm$prior$tau.sq$max <- parm$prior$tau$max^2
	parm$prior$tau.sq$min <- parm$prior$tau$min^2
	parm$prior$inv.tau.sq$max <- 1/parm$prior$tau.sq$min
	parm$prior$inv.tau.sq$min <- 1/parm$prior$tau.sq$max

	parm
	}



########################################

fn.gen.tau  <- function(data, parm)
	{
	###################
	# update tau
	###################

	# only covariates assigned to non-zero row and non-zero group clusters matter

	sum.resid.sq <- 0
	count <- 0

	for (g in 1:parm$clust$G)
		{flag.v <- parm$clust$c.v == g
		z.g.v <- parm$clust$s.mt[,g] > 0

		if ((sum(z.g.v) > 0) & (sum(flag.v)>0))
			{X.g.mt <- parm$X[z.g.v,flag.v]
			a.g.v <- parm$clust$A.mt[z.g.v,g]
			resid.g.mt <- X.g.mt - a.g.v
			sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
			count <- count + sum(z.g.v)*sum(flag.v)
			}

		}

	shape <- parm$prior$tau$alpha + count/2
	rate <- parm$prior$tau$beta + sum.resid.sq/2

	u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
	u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)

      parm$tau <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

	# overwrite to avoid zeros and Inf
	if (round(u.min, digits = 5) == 1) # really close to 1
		{parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$min)
		}
	if (round(u.max, digits = 5) == 0) # really close to 0
		{parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$max)
		}

	###################
	# update tau_int
	###################

	# only covariates assigned to intercept cluster matter

	sum.resid.sq <- 0
	flag.v <- parm$clust$c.v == 0
	count <- parm$clust$C.m0*parm$n2

	if (parm$clust$C.m0>0)
			{X.g.mt <- parm$X[,flag.v]
			a.g.v <- rep(1,parm$n2)
			resid.g.mt <- X.g.mt - a.g.v
			sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
			}


	shape <- parm$prior$tau$alpha + count/2
	rate <- parm$prior$tau$beta + sum.resid.sq/2

	# shares same support as parm$tau
	u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
	u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)

      parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

      parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

	if (gen.u < 1e-5)
		{parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$max)
		}
	if ((1-gen.u) < 1e-5)
		{parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$min)
		}

	parm
	}



fn.gen.tau_0  <- function(data, parm)
	{
	###################
	# update tau_0
	###################

	sum.resid.sq <- 0
	count <- 0

	for (g in 1:parm$clust$G)
		{flag.v <- parm$clust$c.v == g
		z.g.v <- parm$clust$s.mt[,g] > 0

		if ((sum(1-z.g.v) > 0) & (sum(flag.v)>0))
			{X.g.mt <- parm$X[!z.g.v,flag.v]
			resid.g.mt <- X.g.mt
			sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
			count <- count + sum(1-z.g.v)*sum(flag.v)
			}
		}

	shape <- 1 + count/2
	rate <- 1 + sum.resid.sq/2

	# minimum possible value of parm$tau_0 = 1.5 * maximum possible value of parm$tau
	# maximum possible value of parm$tau_0 = 3 * sd(as.vector(data$X))
	u.min <- pgamma(1/9 / var(as.vector(data$X),na.rm=TRUE),shape=shape, rate=rate)
	u.max <- pgamma(1/1.5^2/parm$prior$tau.sq$min,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)

      parm$tau_0 <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

	parm
	}


fn.hyperparameters <- function(data, parm)
	{

	# also updates update tau_int
	parm <- fn.gen.tau(data, parm)

	parm <- fn.gen.tau_0(data, parm)

	parm

	}


fn.funky <- function(s,t)
	{# on log scale
	lgamma(s+t) - lgamma(s)
	}

fn.d <- function(d, parm)
	{

	# formula in review paper by Lijoi and Prunster
	log.lik <- sum(log(parm$b1 + (1:(parm$clust$G-1))*d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-d), (parm$clust$C.m.vec-1)))

	log.lik
	}




fn.poissonDP.hyperparm <- function(data, parm, w=.01, max.d)
	{


	## update parm$d conditional on parm$b1
	## 1/w must be an integer

	d.v <- seq(0,max.d,by=w)
	len <- length(d.v)
	d.v <- d.v[-len]
	len <- len-1

	log.lik.v <- sapply(d.v, fn.d, parm)
	# putting 1/2 prior mass on 0 and remaining spread uniformly on positive points in d.v
	log.p.v <- log(.5) + c(0,  rep(-log(len-1),(len-1)))

	log.post.v <- log.lik.v + log.p.v
	log.post.v <- log.post.v - max(log.post.v)
	post.v <- exp(log.post.v)
	post.v <- post.v/sum(post.v)

	# plot(d.v, post.v, type="l")

	prop.d <- sample(d.v, size=1, prob=post.v)

	if (prop.d > 0)
		{prop.d <- runif(n=1, min=(prop.d-w), max=(prop.d+w))
		}

	if (prop.d != parm$d)
		{
		# MH ratio for independent proposals and
		# prior same for all d (which is true if 0 wp .5 and \in (0,max.d) wp .5)

		log.ratio <- fn.d(d=prop.d, parm) - fn.d(d=parm$d, parm)
		prob <- min(1, exp(log.ratio))
		flip <- rbinom(n=1, size=1, prob=prob)
		if (flip==1)
			{parm$d <- prop.d
			}
		}

	parm

	}

########################################

fn.iter <- function(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm,
                    computeMode)
	{
	parm <- fast_PDP_fn.main(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size, computeMode)

	parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes, computeMode)

	parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)

	parm <- fn.hyperparameters(data, parm)

	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
	parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC: err=",err))
		}

	parm


	}



fn.mcmc <- function(text, true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, dahl=FALSE,
                    computeMode = "R")
	{

	# initialize
	parm <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, computeMode)
	init.parm <- parm

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC at fn.init: err=",err))
		}

	for (cc in 1:n.burn)
		{parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, computeMode)

		if (cc %% 10 == 0)
			{print(paste(text, "BURN = ",cc,date(),"***********"))
			}
		}

	##########################################
	## first get an estimated G cluster
	##########################################


	All.Stuff <- NULL
	#
	All.Stuff$d.v <- All.Stuff$tau_0.v <- All.Stuff$tau.v <- All.Stuff$tau_int.v <- All.Stuff$G.v <- All.Stuff$K.v <- array(,n.reps)
	All.Stuff$row.flip.v  <- array(0,n.reps)
	All.Stuff$nbhd_max <- All.Stuff$col_new_clust.v  <- All.Stuff$col_exit.v <- All.Stuff$col_flip.v  <- array(0,n.reps)

	All.Stuff$pi.mt <- array(0,c(parm$p,parm$p))
	All.Stuff$mean.taxicab.v  <- array(0,n.reps)

	for (cc in 1:n.reps)
		{parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, computeMode)

		All.Stuff$G.v[cc] <- parm$clust$G
		All.Stuff$K.v[cc] <- parm$clust$K
		All.Stuff$tau.v[cc] <- parm$tau
		All.Stuff$tau_0.v[cc] <- parm$tau_0
		All.Stuff$tau_int.v[cc] <- parm$tau_int

		All.Stuff$d.v[cc] <- parm$d

		# summarizing elementwise DP in "fn.groupwise.updates"

		All.Stuff$row.flip.v[cc]  <- parm$clust$row.flip

		All.Stuff$col_new_clust.v[cc]  <- parm$clust$col.new.flag
		All.Stuff$col_flip.v[cc]  <- parm$clust$col.mh.flip
		All.Stuff$col_exit.v[cc]  <- parm$clust$col.mh.exit

		tmp.mat <- array(0,c(parm$p,parm$p))

		for (jj in 1:parm$clust$G)
			{indx.jj <- which(parm$clust$c.v==jj)
			tmp.mat[indx.jj,indx.jj] <- 1
		}

		All.Stuff$pi.mt <- All.Stuff$pi.mt + tmp.mat

		All.Stuff$mean.taxicab.v[cc] <- mean(true_parm$clust$nbhd.matrix != tmp.mat)

		All.Stuff$nbhd_max[cc] <- round(parm$clust$nbhd_max_dist, digits = 2)

		if (cc %% 10 == 0)
			{print(paste(text, "REPS = ",cc,date(),"***********"))
			}

		} # end for loop in cc


	All.Stuff$pi.mt <- All.Stuff$pi.mt/n.reps

	All.Stuff$parm <- parm
	All.Stuff$init.parm <- init.parm

	###

	All.Stuff
	}

