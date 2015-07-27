
PDP_fn.compact.consistency.check <- function(parm)
	{err <- 0

	if (max(sort(unique(parm$clust$c.v))) != parm$clust$G)
		{err <- 1
		}

	if (sum(parm$clust$C.m.vec==0) > 0)
		{err <- 1.5
		}

	err <- PDP_fn.consistency.check(parm)


	err
	}

PDP_fn.check.nbhd <- function(parm)
{
  max_dist.v <- array(,length(parm$clust$col.nbhd.k))

  for (zz in 1:length(parm$clust$col.nbhd.k))
  {
    v1 <- parm$clust$col_post.prob.mt[,parm$clust$col.nbhd.k[zz]]
    m1 <- matrix(parm$clust$col_post.prob.mt[,parm$clust$col.nbhd[[zz]]],ncol=length(parm$clust$col.nbhd[[zz]]))
    max_dist.v[zz] <- max(2*(1-colSums(sqrt(v1*m1))))
  }

  parm$clust$nbhd_max_dist <- max(max_dist.v)

  parm

}

PDP_fn.consistency.check <- function(parm)
	{err <- 0


	if (sum(parm$clust$C.m.vec) + parm$clust$C.m0 != parm$p)
		{err <- 4
		}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 9
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 10
		}


	err
	}


PDP_fn.swap.clusters <- function(parm, g1, g2)
	{

	####################################################
	# swap the group labels g1 and g2
	####################################################

	ind1 <- parm$clust$c.v == g1
	ind2 <- parm$clust$c.v == g2
	parm$clust$c.v[ind1] <- g2
	parm$clust$c.v[ind2] <- g1

	buffer <- parm$clust$s.mt[,g1]
	parm$clust$s.mt[,g1] <- parm$clust$s.mt[,g2] # HOT 3
	parm$clust$s.mt[,g2] <- buffer

	buffer <- parm$clust$beta.v[g1]
	parm$clust$beta.v[g1] <- parm$clust$beta.v[g2]
	parm$clust$beta.v[g2] <- buffer

	buffer <- parm$clust$gamma.v[g1]
	parm$clust$gamma.v[g1] <- parm$clust$gamma.v[g2]
	parm$clust$gamma.v[g2] <- buffer

	buffer <- parm$clust$C.m.vec[g1]
      parm$clust$C.m.vec[g1] <- parm$clust$C.m.vec[g2]
      parm$clust$C.m.vec[g2] <- buffer

	buffer <- parm$clust$small.indx[g1]
      parm$clust$small.indx[g1] <- parm$clust$small.indx[g2]
      parm$clust$small.indx[g2] <- buffer

	buffer <- parm$clust$order.v[g1]
      parm$clust$order.v[g1] <- parm$clust$order.v[g2]
      parm$clust$order.v[g2] <- buffer

	#####################

	buffer <- parm$clust$A.mt[,g1]
	parm$clust$A.mt[,g1] <- parm$clust$A.mt[,g2]
	parm$clust$A.mt[,g2] <- buffer

	buffer <- parm$clust$B.mt[,(g1+1)]
	parm$clust$B.mt[,(g1+1)] <- parm$clust$B.mt[,(g2+1)]
	parm$clust$B.mt[,(g2+1)] <- buffer

	# first swap columns
	buffer <- parm$clust$tBB.mt[,(g1+1)]
	parm$clust$tBB.mt[,(g1+1)] <- parm$clust$tBB.mt[,(g2+1)] # HOT 3
	parm$clust$tBB.mt[,(g2+1)] <- buffer
	# then swap rows
	buffer <- parm$clust$tBB.mt[(g1+1),]
	parm$clust$tBB.mt[(g1+1),] <- parm$clust$tBB.mt[(g2+1),]
	parm$clust$tBB.mt[(g2+1),] <- buffer

	parm
	}

PDP_fn.clip.clusters <- function(parm, keep)
	{

	parm$clust$s.mt <- parm$clust$s.mt[,keep]

	parm$clust$beta.v <- parm$clust$beta.v[keep]
	parm$clust$gamma.v <- parm$clust$gamma.v[keep]
	parm$clust$C.m.vec <- parm$clust$C.m.vec[keep]
      parm$clust$small.indx <- parm$clust$small.indx[keep]
      parm$clust$order.v <- parm$clust$order.v[keep]

     	parm$clust$A.mt <- parm$clust$A.mt[,keep]

	indx2 <- c(1,(keep+1))
	parm$clust$B.mt <- parm$clust$B.mt[,indx2]
	parm$clust$tBB.mt <- parm$clust$tBB.mt[indx2,indx2]

	parm
	}



###########################################################


PDP_fn.log.lik <- function(gg, x.mt, parm, colSums)
	{
	if (gg > 0)
		{a2.v <- parm$clust$A.mt[,gg]
		z.g.v <- parm$clust$s.mt[,gg] > 0
		log.lik.v <- rep(0,ncol(x.mt))

		if (sum(z.g.v) > 0)
			{a2.1.v <- parm$clust$A.mt[z.g.v,gg]
			small.X.1 <- matrix(x.mt[z.g.v,], ncol=ncol(x.mt))
			if (colSums)
				{tmp <- colSums(-.5*(small.X.1 - a2.1.v)^2)
				}
			if (!colSums)
				{tmp <- sum(-.5*(small.X.1 - a2.1.v)^2)
				}
			log.lik.v <- log.lik.v + tmp/parm$tau^2 + sum(z.g.v)*(-.5*log(2*pi)-log(parm$tau))
			}
		if (sum(1-z.g.v) > 0)
			{small.X.0 <- matrix(x.mt[!z.g.v,], ncol=ncol(x.mt))
			if (colSums)
				{tmp <- colSums(-.5*small.X.0^2)
				}
			if (!colSums)
				{tmp <- sum(-.5*small.X.0^2)
				}
			log.lik.v <- log.lik.v + tmp/parm$tau_0^2 + sum(1-z.g.v)*(-.5*log(2*pi)-log(parm$tau_0))
			}
		}

	if (gg == 0)
		{a2.v <- rep(1,parm$n2)
		small.X <- x.mt
		if (colSums)
			{tmp <- colSums(-.5*(small.X - a2.v)^2)
			}
		if (!colSums)
			{tmp <- sum(-.5*(small.X - a2.v)^2)
			}

		log.lik.v <-  tmp/parm$tau_int^2 -parm$n2*.5*log(2*pi)-parm$n2*log(parm$tau_int)
		}

	log.lik.v
	}



PDP_fn.nbhd <- function(relative_I, parm, max.col.nbhd.size)
{
  if (length(relative_I)>1)
  {relative_k <- sample(relative_I, size=1)
  }
  if (length(relative_I)==1)
  {relative_k <- relative_I
  }

  post.prob.mt <- parm$clust$col_subset_post.prob.mt

  tmp1.mt <- matrix(post.prob.mt[,relative_I], ncol=length(relative_I))
  tmp2.v <- post.prob.mt[,relative_k]
  tmp3.mt <- sqrt(tmp1.mt * tmp2.v)
  H.v <-  2*(1-colSums(tmp3.mt))

  cutoff <- parm$col.delta
  flag.v <- which(H.v <= cutoff)
  relative_I.k <- relative_I[flag.v]

  if (length(relative_I.k) > max.col.nbhd.size)
  {cutoff <- quantile(H.v[flag.v], probs=max.col.nbhd.size/length(relative_I.k))
   relative_I.k <- relative_I[which(H.v <= cutoff)]
  }

  relative_I.k <- sort(relative_I.k)

  relative_I <- sort(setdiff(relative_I, relative_I.k))
  relative_I <- sort(relative_I)

  list(relative_k, relative_I.k, relative_I)

}


PDP_fn.post.prob.and.delta <- function(parm, max.col.nbhd.size, computeMode)

{

  col.subset <- 1:parm$p

  if (computeMode$useR) {

    ################################################
    ### Compute pmf of cluster variables w_1,...,w_p
    ###############################################

    prior.prob.v <- c(parm$clust$C.m0, parm$clust$C.m.vec)
    small <- 1e-3 # compared to 1
    prior.prob.v[prior.prob.v < small] <- small

    subset_log.ss.mt <- array(,c((parm$clust$G+1), length(col.subset)))

    for (gg in 0:parm$clust$G)
    {subset_log.ss.mt[(gg+1),] <- PDP_fn.log.lik(gg, x.mt=parm$X[,col.subset], parm, colSums=TRUE)
    }

    subset_log.ss.mt <- subset_log.ss.mt + log(prior.prob.v)

    maxx.v <- apply(subset_log.ss.mt, 2, max)
    subset_log.ss.mt <- t(t(subset_log.ss.mt) - maxx.v)
    subset_ss.mt <- exp(subset_log.ss.mt)

    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    # TODO Ask GS if second normalization is really necessary

    # replace zeros by "small"
    small2 <- 1e-5
    subset_ss.mt[subset_ss.mt < small2] <- small2

    # again normalize
    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    parm$clust$col_post.prob.mt <- array(,c((parm$clust$G+1), parm$p))
    parm$clust$col_post.prob.mt[,col.subset] <- subset_ss.mt

    dimnames(parm$clust$col_post.prob.mt) <- list(0:parm$clust$G, 1:parm$p)

    parm$clust$col_subset_post.prob.mt <- subset_ss.mt
    dimnames(parm$clust$col_subset_post.prob.mt) <- list(0:parm$clust$G, 1:length(col.subset))


    #########################################
    ### now compute the delta-neighborhoods
    #########################################

    # savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only

    parm$clust$col.nbhd <- NULL
    parm$clust$col.nbhd.k <- NULL
    relative_I <- 1:length(col.subset)

    while (length(relative_I)>=1)
    {tmp <- PDP_fn.nbhd(relative_I, parm, max.col.nbhd.size)
     relative_k <- tmp[[1]]
     relative_I.k <- tmp[[2]]
     relative_I <- tmp[[3]]
     #
     parm$clust$col.nbhd <- c(parm$clust$col.nbhd, list(col.subset[relative_I.k]))
     parm$clust$col.nbhd.k <- c(parm$clust$col.nbhd.k, col.subset[relative_k])
    }

  } else { # computeMode != "R"

  ############################################################################
  ### SG did not copy MS's code copied from "element_fn.post.prob.and.delta" here
  ############################################################################

  } # computeMode

  ## END

  ### SG: how good are these nbhds?

  parm <- PDP_fn.check.nbhd(parm)

  parm
}


###########################################################

PDP_fn.gibbs <- function(k, parm, data, computeMode)
{	k <- parm$k

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("GIBBS - 0: failed consistency check: err=",err))
		}

	old.c.k <- parm$clust$c.v[k]

	###############

	if (old.c.k > 0)
		{parm$clust$C.m.vec[old.c.k] <- parm$clust$C.m.vec[old.c.k] - 1
		}
	if (old.c.k == 0)
		{parm$clust$C.m0 <- parm$clust$C.m0 - 1
		}

	x.mt <- matrix(parm$X[,k], ncol=1)

	if (computeMode$useR | computeMode$exactBitStream) {

	# intercept cluster or any existing cluster
	L.v <- sapply(0:parm$clust$G, PDP_fn.log.lik, x.mt, parm, colSums=FALSE) # HOT


  } else { # computeMode
    # engine <- createEngine(sort = TRUE)
    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)
    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

# 	abs <- abs(L.v - test$logLikelihood)
# 	pass <- all(abs < 1e-14)
# 	if(!pass) {
# 	  stop("bad compute")
# 	}

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling
	#######################################################

	emptied.indx <- which(parm$clust$C.m.vec==0)
	new.G <- parm$clust$G - length(emptied.indx)

  if (length(emptied.indx) >0)
  	{
    	new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3

    	if (computeMode$useR) {
    	  new.n.vec <- array(,parm$clust$K)
    	  for (pp in 1:parm$clust$K) {
    	    new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
    	  }
    	} else { # computeMode
    	  new.n.vec <- .fastTabulate(new.s.mt, parm$clust$K)
    	} # computeMode

    	emptied.s.indx <- which(new.n.vec==0)
    	new.K <- parm$clust$K-length(emptied.s.indx)
   	 new.n0 <- sum(new.s.mt==0) # HOT 3
   	}

  if (length(emptied.indx) ==0)
  	{
    	new.s.mt <- parm$clust$s.mt
    	new.n.vec <- parm$clust$n.vec
    	emptied.s.indx <- which(new.n.vec==0)
    	new.K <- parm$clust$K
    	new.n0 <- parm$clust$n0
    }

  ## generate auxilliary P vector

  tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  tmp.alpha <- tmp.M+new.n.vec
  tmp.alpha[emptied.s.indx] <- 0
  P.aux <- rgamma(parm$clust$K+1,c(parm$clust$M0+new.n0,tmp.alpha),1)
  P.aux <- P.aux/sum(P.aux)

  ## marginal likelihood of new cluster
  marg.log.lik.v <- array(,length(x.mt))
  for (tt in 1:length(x.mt))
	{
    	tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) # HOT 3
    	tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
   	marg.log.lik.v[tt] <- log(sum(tmp.lik.v*P.aux))
  	}
  marg.log.lik <- sum(marg.log.lik.v)

  L.v <- c(L.v, marg.log.lik)

	log.prior.v <- array(NA, (2+parm$clust$G))

  # allow -Inf's in Gibbs sampling log-prior (just emptied clusters)
	if (length(emptied.indx) >0)
		{log.prior.v[-(emptied.indx+1)] <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), (parm$b1+new.G*parm$d)))
		log.prior.v[emptied.indx+1] <- -Inf
		}

	if (length(emptied.indx) ==0)
		{log.prior.v <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
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

	new.c.k <- sample(0:(parm$clust$G+1), size=1, replace=TRUE, prob=parm$clust$post.k)
	parm$clust$c.v[k] <- new.c.k
	new.flag <- new.c.k == (parm$clust$G+1)

	#######################

	count.0 <- sum(new.c.k==0)
	parm$clust$C.m0 <- parm$clust$C.m0 + count.0

	if ((count.0 == 0)&(!new.flag))
		{parm$clust$C.m.vec[new.c.k] <- parm$clust$C.m.vec[new.c.k] + 1
		}

	if (new.flag)
  {
    ###generate the latent vector first, condition on the single kth column
    cand.s.v.k <- array(,length(x.mt))
    for (tt in 1:length(x.mt))
	{
      	tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)
      	tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
      	tmp.prob.v <- tmp.lik.v*P.aux
      	prob.gen.v <- tmp.prob.v/sum(tmp.prob.v)
      	cand.s.v.k[tt]<-sample(0:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v)
        }
    parm$cand$s.v.k <- cand.s.v.k

    parm$cand$n.vec.k <- array(,parm$clust$K)
    for (gg in 1:parm$clust$K)
    	{parm$cand$n.vec.k[gg] <- sum(cand.s.v.k==gg)
    	}
    parm$cand$n0.k <- sum(cand.s.v.k==0)

  ##################
   parm$clust$G <- parm$clust$G + 1

		parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)

		parm$clust$beta.v <- c(parm$clust$beta.v, 0)
		parm$clust$gamma.v <- c(parm$clust$gamma.v,0)

   parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)

		parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)

   parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
   parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k

		parm$N <- sum(parm$clust$n.vec) + parm$clust$n0

   tmp.a.v <- array(,parm$n2)
   s.G.v <- parm$cand$s.v.k
		indxx <- s.G.v==0
		tmp.a.v[indxx] <- 0
		tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]]
		#
		parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
		parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)

		if (computeMode$useR) {
		  parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
		} else {
		  parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
		}
  } # end  if (new.flag)


	list(parm, new.flag)
}

###########################################################

PDP_fn.fast_col <- function(cc, parm, data, computeMode)
{
  k <- parm$k <- parm$clust$col.nbhd.k[[cc]]
  I.k <- parm$clust$col.nbhd[[cc]]

   err <- PDP_fn.consistency.check(parm)
   if (err > 0)
   {stop(paste("FAST - 0: failed consistency check: err=",err))
   }

  # store so that we can revert to this state if MH propsal is rejected
  init.cc.parm <- parm

  subset.c <- parm$clust$c.v[I.k]

  # Note to MS: I've ignored the branching conditions
  #             in "elementwise_DP.functions.R" based on computeMode$useR

  parm$clust$C.m.vec.k <- array(,parm$clust$G)
  for (gg in 1:parm$clust$G) {
    parm$clust$C.m.vec.k[gg] <- sum(subset.c==gg)
  }
  parm$clust$C.m0.k <- sum(subset.c==0)

  parm$clust$C.m.vec.k.comp <- parm$clust$C.m.vec - parm$clust$C.m.vec.k
  parm$clust$C.m0.k.comp <- parm$clust$C.m0 - parm$clust$C.m0.k

  old.c.k <- parm$clust$c.v[k]

  x.mt <- matrix(parm$X[,k], ncol=1)

  if (computeMode$useR | computeMode$exactBitStream) {

    # intercept cluster or any existing cluster
    L.v <- sapply(0:parm$clust$G, PDP_fn.log.lik, x.mt, parm, colSums=FALSE) # HOT


  } else { # computeMode
    # engine <- createEngine(sort = TRUE)
    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)
    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

  #   abs <- abs(L.v - test$logLikelihood)
  # 	pass <- all(abs < 1e-14)
  # 	if(!pass) {
  # 	  stop("bad compute")
  # 	}

   #######################################################
   ### emptied clusters are gone forever under Gibbs sampling
   #######################################################

   emptied.indx <- which(parm$clust$C.m.vec.comp==0)
   new.G <- parm$clust$G - length(emptied.indx)

  if (length(emptied.indx) >0)
  {
    new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3

    if (computeMode$useR) {
      new.n.vec <- array(,parm$clust$K)
      for (pp in 1:parm$clust$K) {
        new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
      }
    } else { # computeMode
      new.n.vec <- .fastTabulate(new.s.mt, parm$clust$K)
    } # computeMode

    emptied.s.indx <- which(new.n.vec==0)
    new.K <- parm$clust$K-length(emptied.s.indx)
    new.n0 <- sum(new.s.mt==0) # HOT 3
  }

  if (length(emptied.indx) ==0)
  {
    new.s.mt <- parm$clust$s.mt
    new.n.vec <- parm$clust$n.vec
    emptied.s.indx <- which(new.n.vec==0)
    new.K <- parm$clust$K
    new.n0 <- parm$clust$n0
  }

  ## generate auxilliary P vector

  tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  tmp.alpha <- tmp.M+new.n.vec
  tmp.alpha[emptied.s.indx] <- 0
  P.aux <- rgamma(parm$clust$K+1,c(parm$clust$M0+new.n0,tmp.alpha),1)
  P.aux <- P.aux/sum(P.aux)

  ## marginal likelihood of new cluster
  marg.log.lik.v <- array(,length(x.mt))
  for (tt in 1:length(x.mt))
  {
    tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) # HOT 3
    tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
    marg.log.lik.v[tt] <- log(sum(tmp.lik.v*P.aux))
  }
  marg.log.lik <- sum(marg.log.lik.v)

  L.v <- c(L.v, marg.log.lik)

  ##
   log.prior.v <- array(NA, (2+parm$clust$G))

  spread.mass <- (parm$b1+new.G*parm$d)/(1+length(emptied.indx))

   if (length(emptied.indx) >0)
   {
    log.prior.v[-(emptied.indx+1)] <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
    log.prior.v[emptied.indx+1] <- spread.mass
   }

   if (length(emptied.indx) ==0)
   {log.prior.v <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec-parm$d), spread.mass))
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

   new.c.k <- sample(0:(parm$clust$G+1), size=length(I.k), replace=TRUE, prob=parm$clust$post.k)
   parm$clust$c.v[I.k] <- new.c.k

  exit <- (sum(new.s.k != old.s.k)==0)
  flip <- 1

   new.flag <- new.c.k == (parm$clust$G+1)

   #######################

   count.0 <- sum(new.c.k==0)
   parm$clust$C.m0 <- parm$clust$C.m0 + count.0

   if ((count.0 == 0)&(!new.flag))
   {parm$clust$C.m.vec[new.c.k] <- parm$clust$C.m.vec[new.c.k] + 1
   }

   if (new.flag)
   {
     ###generate the latent vector first, condition on the single kth column
     cand.s.v.k <- array(,length(x.mt))
     for (tt in 1:length(x.mt))
     {
       tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)
       tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
       tmp.prob.v <- tmp.lik.v*P.aux
       prob.gen.v <- tmp.prob.v/sum(tmp.prob.v)
       cand.s.v.k[tt]<-sample(0:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v)
     }
     parm$cand$s.v.k <- cand.s.v.k

     parm$cand$n.vec.k <- array(,parm$clust$K)
     for (gg in 1:parm$clust$K)
     {parm$cand$n.vec.k[gg] <- sum(cand.s.v.k==gg)
     }
     parm$cand$n0.k <- sum(cand.s.v.k==0)

     ##################
     parm$clust$G <- parm$clust$G + 1

     parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)

     parm$clust$beta.v <- c(parm$clust$beta.v, 0)
     parm$clust$gamma.v <- c(parm$clust$gamma.v,0)

     parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)

     parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)

     parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
     parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k

     parm$N <- sum(parm$clust$n.vec) + parm$clust$n0

     tmp.a.v <- array(,parm$n2)
     s.G.v <- parm$cand$s.v.k
     indxx <- s.G.v==0
     tmp.a.v[indxx] <- 0
     tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]]
     #
     parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
     parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)

     if (computeMode$useR) {
       parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
     } else {
       parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
     }
   } # end  if (new.flag)


   list(parm, new.flag)
}


#####################################


PDP_fn.drop <- function(parm)
 	{

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

	parm$clust$G <- sum(parm$clust$C.m.vec>0)
	num.dropped <- sum(parm$clust$C.m.vec==0)

  if (parm$clust$G > 0)
  {
	if (num.dropped > 0)
	{
	for (rr in 1:num.dropped)
		{
		old.label <-  min(which(parm$clust$C.m.vec==0))
		new.label <- max(which(parm$clust$C.m.vec>0))
		stopp <-  max(which(parm$clust$C.m.vec>0)) == parm$clust$G
		if (stopp)
			{break
			}
		parm <- PDP_fn.swap.clusters(parm, g1=new.label, g2=old.label)
 		}
	}

	##########

	keep <- 1:parm$clust$G
	parm <- PDP_fn.clip.clusters(parm, keep)

	###########

	# parm$clust$K does not change (possibly some empty elementwise clusters)

	parm$N <- parm$n2 * parm$clust$G
	parm$clust$s.v <- as.vector(parm$clust$s.mt)

	parm$clust$n0 <- sum(parm$clust$s.v==0)
	parm$clust$n.vec <- array(,parm$clust$K)
	for (ss in 1:parm$clust$K)
		{parm$clust$n.vec[ss] <- sum(parm$clust$s.v==ss)
		}

  }
	parm
	}





fast_PDP_fn.main <- function(parm, data, col.frac.probes, max.col.nbhd.size, computeMode)
{
  p <- parm$p

	##########################
	# compute delta-neighborhoods
	#########################

  # SG: Notice that col.frac.probes is probability of updating the
  #     the delta neighborhoods. It is also used to pick the delta-nbhds to
  #     update

  col_flip <- as.logical(rbinom(n=1,size=1,prob=col.frac.probes))
  if (is.null(parm$clust$col.nbhd.k) | col_flip){
    parm <- PDP_fn.post.prob.and.delta(parm, max.col.nbhd.size, computeMode)
  }


	if (col.frac.probes < 1)
	{parm$subset_nbhd.indx <- sort(sample(1:length(parm$clust$col.nbhd.k), size=round(col.frac.probes*length(parm$clust$col.nbhd.k))))
	}

	if (col.frac.probes == 1)
	{parm$subset_nbhd.indx <- 1:length(parm$clust$col.nbhd.k)
	}


	new.flag.v <- NULL

     for (cc in parm$subset_nbhd.indx)
			{previous.parm <- parm

			parm$k <- parm$clust$col.nbhd.k[[cc]]

      if (length(parm$clust$col.nbhd[[cc]])==1)
      { tmp <- PDP_fn.gibbs(k=parm$k, parm, data, computeMode)
        parm <- tmp[[1]]
        new.flag.v <- c(new.flag.v, tmp[[2]])
      }

			if (length(parm$clust$col.nbhd[[cc]])>1)
			{
			  #tmp <- PDP_fn.fast_col(cc, parm, data, computeMode)
			  #parm <- tmp[[1]]
			  #new.flag.v <- c(new.flag.v, tmp[[2]])
			}

			}

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
			{stop(paste("LOOP: failed consistency check: err=",err))
			}

	parm$clust$col.new.flag <- mean(new.flag.v)

	##########################################
	## Drop empty group clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

	parm <- PDP_fn.drop(parm)

	# now drop empty elementwise clusters

	parm <- element_fn.drop(parm)

	err <- PDP_fn.compact.consistency.check(parm)
	if (err > 0)
		{stop(paste("MAIN FUNCTION END: failed consistency check: err=",err))
		}

	parm
	}



