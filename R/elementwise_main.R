
fn1.update.element.objects <- function(parm, computeMode)
	{

	parm$clust$s.mt <- array(parm$clust$s.v, c(parm$n2, parm$clust$G))

	for (g in 1:parm$clust$G)
		{s.g.v <- parm$clust$s.mt[,g]
		s.pos.indx <- s.g.v > 0
		#
		if (sum(s.pos.indx) > 0)
			{parm$clust$A.mt[s.pos.indx,g] <- parm$clust$phi.v[s.g.v[s.pos.indx]]
			}
		if ((parm$n2-sum(s.pos.indx)) > 0)
			{parm$clust$A.mt[!s.pos.indx,g] <- 0
			}
		parm$clust$B.mt[,(g+1)] <- parm$clust$A.mt[,g]
		}

	parm$clust$theta.v <- as.vector(parm$clust$A.mt)

	if (computeMode$computeR) {

	  parm$clust$n.vec <- array(,parm$clust$K)
	  for (s in 1:parm$clust$K)
		  {parm$clust$n.vec[s] <- sum(parm$clust$s.v == s)
	  }
	  parm$clust$n0 <- sum(parm$clust$s.v == 0)
	}

	if (computeMode$computeC) {

	  all.n.vec <- .fastTabulateVector(parm$clust$s.v, parm$clust$K, TRUE)

	  if (computeMode$computeR) {
      assertEqual(parm$clust$n0, all.n.vec[1])
	    assertEqual(parm$clust$n.vec, all.n.vec[-1])
	  }

	  parm$clust$n0 <- all.n.vec[1]
	  parm$clust$n.vec <- all.n.vec[-1]
	}

	parm
	}


fn2.update.element.objects <- function(parm, computeMode)
	{

	parm$Y <- parm$X.sd <- array(,c(parm$n2, parm$clust$G))

	# group covariate tells which parm$clust$rho.g
	# to use for likelihood calculation
	parm$g <- rep(1:parm$clust$G,each = parm$n2)

	for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 m.g <- parm$clust$C.m.vec[g]

		x.g.v <- x.tmp <- parm$X[,I.g]
		x2.g.v <- x.g.v^2
		if (m.g > 1)
			{x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
			}
		 parm$Y[,g] <- x.g.v

		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1))
			}
		parm$X.sd[,g] <- sd.g.v
		}

	parm$N <- parm$n2 * parm$clust$G

	parm$Y <- as.vector(parm$Y)

	parm$X.sd <- as.vector(parm$X.sd)

	#####################

	parm <- fn1.update.element.objects(parm, computeMode)

	parm
	}



fn.element.DP <- function(data, parm, max.row.nbhd.size, row.frac.probes,
                          computeMode)
{

  if (parm$standardize.X)
  {parm <- fn.standardize_orient.X(parm)}

	# essentially, a Bush-Mac move: given groups, the parm$N=n2XG number of
	# invidividual elements (summaries of microarray elements) belonging to group g>0
	# are updated for s (phi) and z

	parm <- fn2.update.element.objects(parm, computeMode)

	parm <- element_fn.fast.DP(parm, max.row.nbhd.size, row.frac.probes, computeMode)

	#############################
	## Important: do not remove call to fn1.update.element.objects
	## updates A.mt, theta.v, B.mt, tBB.mt, s.v, s.mt, n.vec, n0
	#############################

	parm <- fn1.update.element.objects(parm, computeMode)

  	parm
}



