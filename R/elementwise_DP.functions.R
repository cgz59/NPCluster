

element_fn.compact_consistency.check <- function(parm)
	{err <- 0

	if (max(unique(parm$clust$s.v)) != parm$clust$K)
		{err <- 1
		}
	
	err <- element_fn.consistency.check(parm)

	err
	}


element_fn.consistency.check <- function(parm)
	{err <- 0
	
	if (!is.null(parm$clust$row.nbhd))
		{if (sum(unlist(lapply(parm$clust$row.nbhd, length))) != length(parm$row.subset.I))
			{err <- 2
			}
		}

	if (sum(parm$clust$n.vec) + parm$clust$n0 != parm$N)
		{err <- 4
		}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 5
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 6
		}

	err
	}


element_fn.swap.clusters <- function(parm, g1, g2)
	{

	####################################################
	# swap the clusters with labels g1 and g2
	####################################################

	ind1 <- parm$clust$s.v == g1
	ind2 <- parm$clust$s.v == g2
	parm$clust$s.v[ind1] <- g2
	parm$clust$s.v[ind2] <- g1

	buffer <- parm$clust$phi.v[g1]
	parm$clust$phi.v[g1] <- parm$clust$phi.v[g2]
	parm$clust$phi.v[g2] <- buffer
		
	buffer <- parm$clust$n.vec[g1]
      parm$clust$n.vec[g1] <- parm$clust$n.vec[g2]
      parm$clust$n.vec[g2] <- buffer

	#####################

	parm
	}




########################################################


element_fn.log.lik <- function(mean, sd, num, Y, X.sd)
	{log.lik <- -num/2*log(2*pi*sd^2) - .5*num*(Y-mean)^2/sd^2 - .5*(num-1)*X.sd^2/sd^2 

	log.lik
	}




element_fn.gen.phi <- function(parm, Y, lik.sd)
{ 	
	if (is.null(Y))
		{lik.sd <- Inf 
		Y <- 0 #dummy
		}


	post.sd  <- 1/sqrt(1/parm$clust$tau2^2 + 1/lik.sd^2)

	post.mean <- post.sd^2 * (parm$clust$mu2/parm$clust$tau2^2 + Y/lik.sd^2)

	new.phi <- rnorm(n=1, mean=post.mean, sd=post.sd)

	new.phi
}



element_fn.nbhd <- function(I, parm, max.row.nbhd.size)
	{if (length(I)>1) 
		{k <- sample(I, size=1)
		}
	if (length(I)==1) 
		{k <- I
		}

       H.v <- rep(Inf,parm$N) 

	 post.prob.mt <- parm$clust$post.prob.mt

	 tmp1.mt <- matrix(post.prob.mt[,I], ncol=length(I))
	 tmp2.v <- post.prob.mt[,k]
	 tmp3.mt <- sqrt(tmp1.mt * tmp2.v)
	 H.v[I] <-  2*(1-colSums(tmp3.mt))

	 cutoff <- parm$row.delta 
	 I.k <- which(H.v <= cutoff)

	if (length(I.k) > max.row.nbhd.size)
		{cutoff <- quantile(H.v[I.k], probs=max.row.nbhd.size/length(I.k))
	 	 I.k <- which(H.v <= cutoff)
		}

	 I.k <- sort(I.k)

       I <- sort(setdiff(I, I.k))
       I <- sort(I)

       list(k, I.k, I)

	}



element_fn.post.prob.and.delta <- function(parm, max.row.nbhd.size)

	{

	################################################
	### Compute pmf of cluster variables w_1,...,w_n
	###############################################

	prior.prob.v <- c(parm$clust$n0, parm$clust$n.vec)
	small <- 1e-3 # compared to 1
	prior.prob.v[prior.prob.v < small] <- small

	log.ss.mt <- array(,c(parm$clust$K, parm$N))

	#
	Y.v <- parm$Y[parm$row.subset.I]
	X.sd.v <- parm$X.sd[parm$row.subset.I]
	## g.v also equals parm$g[parm$row.subset.I]  
	g.v <- (parm$row.subset.I-1) %/% parm$n2 + 1
	num.k.v <-  parm$clust$C.m.vec[g.v]

	for (ss in 1:parm$clust$K)
		{log.ss.mt[ss,parm$row.subset.I] <- element_fn.log.lik(mean=parm$clust$phi.v[ss], sd=parm$tau, num=num.k.v, Y=Y.v, X.sd=X.sd.v)
		}

	## adding the row on top corresponding to s=0
	tmp.v <- array(,parm$N)
	tmp.v[parm$row.subset.I] <- element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k.v, Y=Y.v, X.sd=X.sd.v)
	log.ss.mt <- rbind(tmp.v, log.ss.mt)

	log.ss.mt <- log.ss.mt + log(prior.prob.v)

	dimnames(log.ss.mt) <- list(0:parm$clust$K, 1:parm$N)

	maxx.v <- apply(log.ss.mt, 2, max)
	log.ss.mt <- t(t(log.ss.mt) - maxx.v)
	ss.mt <- exp(log.ss.mt)

	col.sums.v <- colSums(ss.mt)
	ss.mt <- t(t(ss.mt)/col.sums.v)

	# replace zeros by "small"
	small <- 1e-5
	ss.mt[ss.mt < small] <- small

	# again normalize 
	col.sums.v <- colSums(ss.mt)
	ss.mt <- t(t(ss.mt)/col.sums.v)

	parm$clust$post.prob.mt <- ss.mt


	#########################################
	### now compute the delta-neighborhoods
	#########################################

	I <- parm$row.subset.I 
	parm$clust$row.nbhd <- NULL
	parm$clust$row.nbhd.k <- NULL
	
	while (length(I)>=1)
		{tmp <- element_fn.nbhd(I, parm, max.row.nbhd.size)
		k <- tmp[[1]]
		I.k <- tmp[[2]]
		I <- tmp[[3]]
		#
		parm$clust$row.nbhd <- c(parm$clust$row.nbhd, list(I.k))
		parm$clust$row.nbhd.k <- c(parm$clust$row.nbhd.k, k)
		}

	parm
	}
	 


###########################################################

element_fn.row.gibbs.DP <- function(parm)
{	
     k <- parm$k
	I.k <- parm$I.k
	s.k <- parm$clust$s.v[k]

	###############
	
	err <- element_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("GIBBS STEP - 0: failed consistency check: err=",err))
		}

	# store just in case
	init.cc.parm <- parm

	###############

	if (s.k > 0)
		{parm$clust$n.vec[s.k] <- parm$clust$n.vec[s.k] - 1
		}

	if (s.k == 0)
		{parm$clust$n0 <- parm$clust$n0 - 1
		}
	
	# sum(parm$clust$n.vec) + parm$clust$n0 == (parm$N-1)

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling move
	#######################################################

	Y.k <- parm$Y[k]
	X.sd.k <- parm$X.sd[k]
	g.k <- (k-1) %/% parm$n2 + 1
	num.k <-  parm$clust$C.m.vec[g.k]

	L.v <- sapply(parm$clust$phi.v, element_fn.log.lik, sd=parm$tau, num=num.k, Y=Y.k, X.sd=X.sd.k) 
	# add possibility that s=0 in front
	L.v <- c(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k, Y=Y.k, X.sd=X.sd.k), L.v)

	new.phi <- element_fn.gen.phi(parm, Y=Y.k, lik.sd=parm$tau/sqrt(num.k))
	# last element is for new cluster
	parm$clust$n.vec <- c(parm$clust$n.vec, 0)
	parm$clust$phi.v <- c(parm$clust$phi.v, new.phi)
	parm$clust$K <- parm$clust$K + 1

	new.L <- element_fn.log.lik(mean=new.phi, sd=parm$tau, num=num.k, Y=Y.k, X.sd=X.sd.k) + dnorm(new.phi, mean=parm$clust$mu2, sd=parm$clust$tau2, log=TRUE)
	L.v <- c(L.v, new.L)

	############
	
	# forcing singletons to leave empty clusters
	log.prior.v <- log(c((parm$clust$M0+parm$clust$n0), parm$clust$n.vec[-parm$clust$K], parm$clust$M))
	
	tmp2 <- log.prior.v + L.v 
	maxx <- max(tmp2)
	tmp2 <- tmp2 - maxx
	tmp2 <- exp(tmp2)

	parm$clust$post.k <- tmp2
	parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)


	########################################################################
	# store current state 
	old.parm <- parm

	########################################################################

	new.s.k <- sample(0:parm$clust$K, size=length(I.k), replace=TRUE, prob=parm$clust$post.k)
				
	parm$clust$s.v[k] <- new.s.k
	if (new.s.k > 0)
		{parm$clust$n.vec[new.s.k] <- parm$clust$n.vec[new.s.k] + 1
		}
	if (new.s.k == 0)
		{parm$clust$n0 <- parm$clust$n0 + 1
		}

	# sum(parm$clust$n.vec) + parm$clust$n0 == parm$N

	err <- element_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("GIBBS STEP - 1: failed consistency check: err=",err))
		}

	parm
}


###########################################################

element_fn.fast.DP.iter <- function(parm)
{	
     k <- parm$k
	I.k <- parm$I.k

	###############
	
	err <- element_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("FORWARD STEP - 0: failed consistency check: err=",err))
		}

	# store so that we can revert to this state if MH propsal is rejected
	init.cc.parm <- parm

	###############

	subset.s <- parm$clust$s.v[I.k]
	parm$clust$n.vec.k <- array(,parm$clust$K)
	for (gg in 1:parm$clust$K)
		{parm$clust$n.vec.k[gg] <- sum(subset.s==gg)
		}
	parm$clust$n0.k <- sum(subset.s==0)

	parm$clust$n.vec.k.comp <- parm$clust$n.vec - parm$clust$n.vec.k
	parm$clust$n0.k.comp <- parm$clust$n0 - parm$clust$n0.k

	#

	# sum(parm$clust$n.vec.k) + parm$clust$n0.k == length(I.k)
	# sum(parm$clust$n.vec.k.comp) + parm$clust$n0.k.comp == (parm$N-length(I.k))

	Y.k <- parm$Y[k]
	X.sd.k <- parm$X.sd[k]
	g.k <- (k-1) %/% parm$n2 + 1
	num.k <-  parm$clust$C.m.vec[g.k]

	L.v <- sapply(parm$clust$phi.v, element_fn.log.lik, sd=parm$tau, num=num.k, Y=Y.k, X.sd=X.sd.k) 
	# add possibility that s=0 in front
	L.v <- c(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k, Y=Y.k, X.sd=X.sd.k), L.v)

	############
	
	log.prior.v <- log(c((parm$clust$M0+parm$clust$n0.k.comp), parm$clust$n.vec.k.comp))

	#######################################################
	# Buffering empty clusters with positive masses using Neal's method
	# NECESSARY under MH algorithm otherwise reverse proposals are impossible
	# and proposal always accepted (bad)
	#######################################################

	newly.empty.indx <- 1+which(parm$clust$n.vec.k.comp==0)
	
	# possibly overwrite some 0's 
	if (length(newly.empty.indx) > 0)
		{log.prior.v[newly.empty.indx] <- log(parm$clust$M/length(newly.empty.indx))
		}

	tmp2 <- log.prior.v + L.v 
	maxx <- max(tmp2)
	tmp2 <- tmp2 - maxx
	tmp2 <- exp(tmp2)

	parm$clust$post.k <- tmp2
	parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)


	########################################################################
	# store current state 
	old.parm <- parm

	########################################################################

	new.s.k <- sample(0:parm$clust$K, size=length(I.k), replace=TRUE, prob=parm$clust$post.k)
	old.s.k <- old.parm$clust$s.v[I.k] 

	exit <- (sum(new.s.k != old.s.k)==0)
	flip <- 1

  if (!exit) # CONTINUE W/O EXITING FUNCTION
  {				
	parm$clust$s.v[I.k] <- new.s.k
	
	rho.prop <- 0
	parm$clust$n.vec <- array(0, parm$clust$K)
	
	for (gg in 0:parm$clust$K)
		{flag.gg <- (new.s.k==gg)
		 count.gg <- sum(flag.gg)

		if (gg > 0)
			{parm$clust$n.vec[gg] <- parm$clust$n.vec.k.comp[gg] + count.gg
			}

		if (gg == 0)
			{parm$clust$n0 <- parm$clust$n0.k.comp + count.gg
			}
			 
		 if (count.gg > 0)
			{rho.prop <- rho.prop + log(parm$clust$post.k[gg+1])*count.gg
			}
		}

	# sum(parm$clust$n.vec) + parm$clust$n0 == parm$N

			
	########################################## 
	########################################## 
	##### Computing proposal prob of reverse move
	########################################## 
	########################################## 

	for (gg in 0:parm$clust$K)
		{flag.gg <- (old.s.k==gg)
		 count.gg <- sum(flag.gg)

		 if (count.gg > 0)
			{rho.prop <- rho.prop - log(old.parm$clust$post.k[gg+1])*count.gg
			}
		}


	#######################################################
	#######################################################
	########## computing true log-ratio
	#######################################################
	#######################################################

	# formula on page 4 of 11/13/11 notes

	new.indx <- which(parm$clust$n.vec>0)
	old.indx <- which(old.parm$clust$n.vec>0)

	new.K <- length(new.indx)
	old.K <- length(old.indx)

	new.is.s0 <- parm$clust$n0 > 0
	old.is.s0 <- old.parm$clust$n0 > 0

	rho.tru <- sum(lgamma(parm$clust$n.vec[new.indx])) - sum(lgamma(old.parm$clust$n.vec[old.indx])) + log(parm$clust$M) * (new.K - old.K)
	rho.tru <- rho.tru + sum(lgamma(parm$clust$M0 + parm$clust$n0)) - sum(lgamma(parm$clust$M0 + old.parm$clust$n0)) + log(parm$clust$M0) * (new.is.s0 - old.is.s0)

	###################

	Y.k.v <- parm$Y[I.k]
	X.sd.k.v <- parm$X.sd[I.k]
	## g.k.v also equals parm$g[I.k]  
	g.k.v <- (I.k-1) %/% parm$n2 + 1
	num.k.v <-  parm$clust$C.m.vec[g.k.v]

	new.log.lik <- old.log.lik <- 0

	new.s.pos.indx <- new.s.k[new.s.k > 0]
	old.s.pos.indx <- old.s.k[old.s.k > 0]

	if (length(new.s.pos.indx) > 0)
		{new.log.lik <- new.log.lik + sum(element_fn.log.lik(mean=parm$clust$phi.v[new.s.pos.indx], sd=parm$tau, num=num.k.v[new.s.k > 0], Y=Y.k.v[new.s.k > 0], X.sd=X.sd.k.v[new.s.k > 0]))
		}
	if (length(old.s.pos.indx) > 0)
		{old.log.lik <- old.log.lik + sum(element_fn.log.lik(mean=old.parm$clust$phi.v[old.s.pos.indx], sd=parm$tau, num=num.k.v[old.s.k > 0], Y=Y.k.v[old.s.k > 0], X.sd=X.sd.k.v[old.s.k > 0]))
		}

	if (sum(new.s.k == 0)>0)
		{new.log.lik <- new.log.lik + sum(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k.v[new.s.k == 0], Y=Y.k.v[new.s.k == 0], X.sd=X.sd.k.v[new.s.k == 0]))
		}
	if (sum(old.s.k == 0)>0)
		{old.log.lik <- old.log.lik + sum(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k.v[old.s.k == 0], Y=Y.k.v[old.s.k == 0], X.sd=X.sd.k.v[old.s.k == 0]))
		}
	
	rho.tru <- rho.tru + new.log.lik - old.log.lik

	########## toss a coin #################
	
	prob <- exp(min((rho.tru-rho.prop),0))
	flip <- as.logical(rbinom(n=1, size=1, prob=prob))
	if (!flip) {parm <- init.cc.parm}

 } # end BIG if (!exit) loop


	list(parm, flip)
}

########################################

element_fn.gen.new <- function(parm, num.gen)
	{
	parm$clust$K.old <- parm$clust$K
	
	for (tt in 1:num.gen)
		{
		# generate from prior
		new.phi <- element_fn.gen.phi(parm, Y=NULL, lik.sd=NULL)
		#	
		parm$clust$n.vec <- c(parm$clust$n.vec, 0)
		parm$clust$phi.v <- c(parm$clust$phi.v, new.phi)
		parm$clust$K <- parm$clust$K + 1
		#
		}

	parm
	}

#####################################


element_fn.drop <- function(parm)
 	{

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$K equal to number of non-empty clusters 
	## (ii) Retain only clusters  1,...,parm$clust$K
	#########################################

	parm$clust$K <- sum(parm$clust$n.vec>0)
	num.dropped <- sum(parm$clust$n.vec==0)

	if (num.dropped > 0)
	{
	for (rr in 1:num.dropped)
		{
		old.label <-  min(which(parm$clust$n.vec==0))
		new.label <- max(which(parm$clust$n.vec>0))
		stopp <-  max(which(parm$clust$n.vec>0)) == parm$clust$K
		if (stopp)
			{break
			}
		parm <- element_fn.swap.clusters(parm, g1=new.label, g2=old.label)
 		}
	} 
		
	##########

	keep <- 1:parm$clust$K

	parm$clust$phi.v <- parm$clust$phi.v[keep]	
	parm$clust$n.vec <- parm$clust$n.vec[keep]
	parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)
	
	parm
	}



element_fn.fast.DP <- function(parm, max.row.nbhd.size, row.frac.probes)
{		

	if (row.frac.probes < 1)
		{parm$row.subset.I <- sort(sample(1:parm$N, round(row.frac.probes*parm$N)))
		}

	if (row.frac.probes == 1)
		{parm$row.subset.I <- 1:parm$N
		}

	##########################
	# compute delta-neighborhoods
	#########################

	parm <- element_fn.post.prob.and.delta(parm, max.row.nbhd.size)

	err <- element_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("ELEMENT_DP MAIN: failed consistency check: err=",err))
		}
	
	#################################
	# generate new empty clusters from prior
	#################################

	parm <- element_fn.gen.new(parm, num.gen=5)

	##########################
	# M-H move for random effects of delta-nieghborhoods
	#########################

	num.nbhds <- length(parm$clust$row.nbhd)
	flip.v <- NULL

      for (cc in 1:num.nbhds)
		{parm$k <- parm$clust$row.nbhd.k[cc]
		parm$I.k <- parm$clust$row.nbhd[[cc]]

		 # parm$I.k will contain at least parm$k

		if (length(parm$I.k) > 1)
			{tmp <- element_fn.fast.DP.iter(parm)
		 	parm <- tmp[[1]]
		 	flip.v <- c(flip.v, tmp[[2]])
			}

		if (length(parm$I.k) == 1)
			{parm <- element_fn.row.gibbs.DP(parm)
			}
		 #
		}


	######################################

	parm$clust$row.flip <- mean(flip.v)

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$K equal to number of non-empty clusters 
	## (ii) Retain only clusters  1,...,parm$clust$K
	#########################################

	parm <- element_fn.drop(parm)

	#################################
	# uprow.data parm$phi.v conditional on clusters
	#################################


	### implements data reduction formula on page 1 of 11/14/11 notes

	indxx <- parm$clust$s.v>0
	num.v <- parm$clust$C.m.vec[parm$g][indxx]
	Y_num.v <- parm$Y[indxx] * num.v
	s.pos <- parm$clust$s.v[indxx]
	Y_num.j.v <- tapply(Y_num.v, s.pos, sum)
	num.j.v <- tapply(num.v,s.pos, sum)
	Y.j.v <- Y_num.j.v/num.j.v


	for (j in 1:parm$clust$K)
		{parm$clust$phi.v[j] <- element_fn.gen.phi(parm, Y=Y.j.v[j], lik.sd=parm$tau/sqrt(num.j.v[j]))
		}

	err <- element_fn.compact_consistency.check(parm)
	if (err > 0)
		{stop(paste("ELEMENT_DP MAIN END: failed consistency check: err=",err))
		}

	parm
	}



