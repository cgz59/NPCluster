

fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
	{

	# Computes multivariate normal density function
	# a little faster than dmvnorm function of R!

	if (missing(inv.sigma))
		{inv.sigma <- solve(sigma)
		}

	logdet <- as.numeric(determinant(inv.sigma, log=TRUE)$mod)
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
		{if (sum(unlist(lapply(parm$clust$col.nbhd, length))) != length(parm$col.subset.I))
			{err <- 1
			}
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

	if ((sum(parm$clust$C.m.vec) + parm$clust$C.m0) != p)
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




fn.init <- function(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes)
	{
	parm <- true_parm

	parm$clust$C.m0 <- 0

	parm$clust$A.mt <- array(,c(n2,parm$clust$G))

	for (g in 1:parm$clust$G)
		{z.v <- parm$clust$s.mt[,g]>0
		parm$clust$A.mt[z.v,g] <- parm$clust$phi.v[parm$clust$s.mt[z.v,g]]
		parm$clust$A.mt[-z.v,g] <- 0
		}

	parm$clust$B.mt <- cbind(rep(1,n2), parm$clust$A.mt)

	parm$shift <- true$shift

	parm$tau_int <- parm$tau

	parm$X <- data$X

	############################
	# For delta neighborhoods
	############################

	parm$col.delta <- .05

	# delta-neighborhood threshold for elements
	parm$row.delta <- .1
	
	parm

	}



########################################

fn.iter <- function(data, parm, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm)
	{
	parm <- PDP_fn.main(parm, data, col.frac.probes)

	parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes)

	parm$clust$B.mt <- cbind(rep(1,n2), parm$clust$A.mt)
	parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC: err=",err))
		}

	parm


	}



fn.mcmc <- function(text, true, data, n.burn, n.reps, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm)
	{

	# initialize
	parm <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes)
	init.parm <- parm

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC at fn.init: err=",err))
		}

	for (cc in 1:n.burn)
		{parm <- fn.iter(data, parm, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm)

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
	All.Stuff$gibbs.new_col_clust.v  <- All.Stuff$col_flip.v  <- array(0,n.reps)

	All.Stuff$mean.taxicab.v  <- array(0,n.reps)

	for (cc in 1:n.reps)
		{parm <- fn.iter(data, parm, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm)

		All.Stuff$G.v[cc] <- parm$clust$G
		All.Stuff$K.v[cc] <- parm$clust$K
		All.Stuff$tau.v[cc] <- parm$tau 
		All.Stuff$tau_0.v[cc] <- parm$tau_0
		All.Stuff$tau_int.v[cc] <- parm$tau_int

		All.Stuff$d.v[cc] <- parm$d

		# summarizing elementwise DP in "fn.groupwise.updates"

		All.Stuff$row.flip.v[cc]  <- parm$clust$row.flip

		All.Stuff$gibbs.new_col_clust.v[cc]  <- parm$clust$gibbs.new.flag
		
		tmp.mat <- array(0,c(p,p))

		for (jj in 1:parm$clust$G)
			{indx.jj <- which(parm$clust$c.v==jj)
			tmp.mat[indx.jj,indx.jj] <- 1
			}

		All.Stuff$mean.taxicab.v[cc] <- mean(true_parm$clust$nbhd.matrix != tmp.mat)

		if (cc %% 10 == 0)
			{print(paste(text, "REPS = ",cc,date(),"***********"))
			}

		} # end for loop in cc


	All.Stuff$parm <- parm
	All.Stuff$init.parm <- init.parm

	###
	
	All.Stuff
	}

