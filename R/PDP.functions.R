



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


PDP_fn.consistency.check <- function(parm)
	{err <- 0
	
	
	if (sum(parm$clust$C.m.vec) + parm$clust$C.m0 != p)
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
	parm$clust$s.mt[,g1] <- parm$clust$s.mt[,g2]
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
	parm$clust$tBB.mt[,(g1+1)] <- parm$clust$tBB.mt[,(g2+1)]
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
		{a2.v <- rep(1,n2)
		small.X <- x.mt
		if (colSums)
			{tmp <- colSums(-.5*(small.X - a2.v)^2)
			}
		if (!colSums)
			{tmp <- sum(-.5*(small.X - a2.v)^2)
			}

		log.lik.v <-  tmp/parm$tau_int^2 -n2*.5*log(2*pi)-n2*log(parm$tau_int)
		}

	log.lik.v
	}



PDP_fn.gen.new.column <- function(I.k, in.parm)
	{new.parm <- in.parm

	## correctly set these objects:
	## new.parm$N, new.parm$clust$n.vec, new.parm$clust$n0

	#############################################
	# copied from fn2.update.row.objects
	#############################################

	m.g <- length(I.k)

	x.g.v <- x.tmp <- new.parm$X[,I.k]
	x2.g.v <- x.g.v^2
	if (m.g > 1)
		{x.g.v <- rowMeans(x.g.v)
		x2.g.v <- rowMeans(x.tmp^2)
		}

	sd.g.v <- rep(0,n2)
	if (m.g > 1)
		{sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1))
		}

	prior.prob.v <- c((new.parm$clust$n0+new.parm$clust$M0), (new.parm$clust$n.vec+new.parm$clust$M/new.parm$clust$K))
	small <- 1e-3 # compared to 1
	prior.prob.v[prior.prob.v < small] <- small

	log.ss.mt <- array(,c(new.parm$clust$K, n2))

	for (ss in 1:new.parm$clust$K)
		{log.ss.mt[ss,] <- element_fn.log.lik(mean=new.parm$clust$phi.v[ss], sd=new.parm$tau, num=length(I.k), Y=x.g.v, X.sd=sd.g.v)
		}

	## adding the row on top corresponding to s=0
	tmp.v <- element_fn.log.lik(mean=0, sd=new.parm$tau_0, num=length(I.k), Y=x.g.v, X.sd=sd.g.v)

	log.ss.mt <- rbind(tmp.v, log.ss.mt)

	log.ss.mt <- log.ss.mt + log(prior.prob.v)

	dimnames(log.ss.mt) <- list(0:new.parm$clust$K, 1:n2)

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

	new.parm$clust$post.prob.mt <- ss.mt
	new.parm$clust$log.post.prob.mt <- log(new.parm$clust$post.prob.mt)

	##############################

	cum.ss.mt <- apply(ss.mt, 2, cumsum)
	dimnames(cum.ss.mt) <- NULL

	u.v <- runif(n=n2)

	tmp.mt <- t(t(cum.ss.mt) > u.v)

	s.v <- new.parm$clust$K + 1 - as.vector(colSums(tmp.mt)) 

	########################

	new.parm$new$s.v.k <- s.v

	new.parm$new$n.vec.k <- array(,new.parm$clust$K)
	for (gg in 1:new.parm$clust$K)
		{new.parm$new$n.vec.k[gg] <- sum(s.v==gg)
		}
	new.parm$new$n0.k <- sum(s.v==0)

	########################

	log.prop.v <- array(,n2)
	for (i in 1:n2)
		{log.prop.v[i] <- new.parm$clust$log.post.prob.mt[(s.v[i]+1),i]
		}

	##############################

	log.lik <- 0
	for (j in 0:new.parm$clust$K)
		{indx.j <- s.v==j
		if ((sum(indx.j)>0)&(j>0))
			{log.lik <- log.lik + sum(dnorm(x.g.v[indx.j], mean=new.parm$clust$phi.v[j],sd=new.parm$tau,log=TRUE))
			}
		if ((sum(indx.j)>0)&(j==0))
			{log.lik <- log.lik + sum(dnorm(x.g.v[indx.j], mean=0,sd=new.parm$tau_0,log=TRUE))
			}
		}

	add.count <- new.parm$clust$M/new.parm$clust$K

	big.joint <- new.parm$clust$K*log(new.parm$clust$M) + as.numeric((new.parm$clust$n0 + new.parm$new$n0.k)>0)*log(new.parm$clust$M0) 
	big.joint <- big.joint + lgamma(new.parm$clust$M0+new.parm$clust$M) - lgamma(new.parm$clust$M0+new.parm$clust$M+n2+new.parm$N) 
	big.joint <- big.joint + lgamma(new.parm$clust$M0+new.parm$clust$n0+new.parm$new$n0.k) - lgamma(new.parm$clust$M0+1) 
	big.joint <- big.joint + sum(lgamma(new.parm$clust$n.vec+new.parm$new$n.vec.k+add.count)) + sum(dnorm(new.parm$clust$phi.v,log=TRUE))

	small.joint <- 0

	small.joint <- new.parm$clust$K*log(new.parm$clust$M) + as.numeric(new.parm$clust$n0>0)*log(new.parm$clust$M0) 
	small.joint <- small.joint + lgamma(new.parm$clust$M0+new.parm$clust$M) - lgamma(new.parm$clust$M0+new.parm$clust$M+new.parm$N) 
	small.joint <- small.joint + lgamma(new.parm$clust$M0+new.parm$clust$n0) - lgamma(new.parm$clust$M0+1) 
	small.joint <- small.joint + sum(lgamma(new.parm$clust$n.vec+add.count)) + sum(dnorm(new.parm$clust$phi.v,log=TRUE))

	log.true <- log.lik + big.joint - small.joint

	new.parm$new$log.w <- log.true - sum(log.prop.v)
	new.parm$new$log.true <- log.true
	new.parm$new$log.prop.v <- log.prop.v

	new.parm$new$log.lik <- log.lik
	new.parm$new$big.joint <- big.joint
	new.parm$new$small.joint <- small.joint

	###########################
		
	new.parm
	}



###########################################################


PDP_fn.gibbs <- function(k, parm, data)
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

	##############################
	
	x.mt <- matrix(parm$X[,k], ncol=1)
	# intercept cluster or any existing cluster
	L.v <- sapply(0:parm$clust$G, PDP_fn.log.lik, x.mt, parm, colSums=FALSE) 

	# create new potential cluster
	# and compute parm$new$log.w
	parm <- PDP_fn.gen.new.column(I.k=k, in.parm=parm)

	L.v <- c(L.v, parm$new$log.w)

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling 
	#######################################################
	
	# allow -Inf's in Gibbs sampling log-prior (just emptied clusters)
	emptied.indx <- which(parm$clust$C.m.vec==0)
	new.G <- parm$clust$G - length(emptied.indx)

	log.prior.v <- array(NA, (2+parm$clust$G))
	
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
		{parm$clust$G <- parm$clust$G + 1
					
		parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)

		parm$clust$beta.v <- c(parm$clust$beta.v, 0)
		parm$clust$gamma.v <- c(parm$clust$gamma.v,0)

		parm$clust$s.v <- c(parm$clust$s.v, parm$new$s.v.k)

		parm$clust$s.mt <- matrix(parm$clust$s.v, nrow=n2)
	
		parm$clust$n.vec <- parm$clust$n.vec + parm$new$n.vec.k
		parm$clust$n0 <- parm$clust$n0 + parm$new$n0.k
	
		parm$N <- sum(parm$clust$n.vec) + parm$clust$n0

		tmp.a.v <- array(,n2)
		s.G.v <- parm$new$s.v.k
		indxx <- s.G.v==0
		tmp.a.v[indxx] <- 0
		tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]]
		#
		parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v)
		parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)
		parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
		}


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

	parm$N <- n2 * parm$clust$G
	parm$clust$s.v <- as.vector(parm$clust$s.mt)

	parm$clust$n0 <- sum(parm$clust$s.v==0)
	parm$clust$n.vec <- array(,parm$clust$K)
	for (ss in 1:parm$clust$K)
		{parm$clust$n.vec[ss] <- sum(parm$clust$s.v==ss)
		}

  }
	parm
	}





PDP_fn.main <- function(parm, data, col.frac.probes)
{		

	if (col.frac.probes < 1)
		{parm$col.subset.I <- sort(sample(1:p, round(col.frac.probes*p)))
		}

	if (col.frac.probes == 1)
		{parm$col.subset.I <- 1:p
		}

	gibbs.new.flag.v <- NULL

     for (cc in parm$col.subset.I)
			{previous.parm <- parm

			parm$k <- cc

			tmp <- PDP_fn.gibbs(k=parm$k, parm, data)
			parm <- tmp[[1]]
			gibbs.new.flag.v <- c(gibbs.new.flag.v, tmp[[2]])
				
			}

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
			{stop(paste("LOOP: failed consistency check: err=",err))
			}		
			
	parm$clust$gibbs.new.flag <- mean(gibbs.new.flag.v)
	
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



