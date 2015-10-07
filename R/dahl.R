
fn.dahl <- function(All.Stuff, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, computeMode = "R")

  {
    parm <- All.Stuff$parm

  All.Stuff$dahl$min <- Inf

  c.matrix <- array(0,c(n.reps,parm$p))

  G.matrix <- posterior$pi.mt

  for (cc in 1:n.reps)
    {
    ### Dahl calculations

    tmp.mat <- array(0,c(p,p))

    for (jj in 0:All.Stuff$G.v[cc])
      {indx.jj <- which(c.matrix[cc,]==jj)
      tmp.mat[indx.jj,indx.jj] <- 1
      }

    dahl.dist <- sum((G.matrix - tmp.mat)^2)

    if (dahl.dist < All.Stuff$dahl$min)
      {All.Stuff$dahl$min <- dahl.dist
    All.Stuff$dahl$G.clust <- c.matrix[cc,]
    All.Stuff$dahl$G <- All.Stuff$G.v[cc]
    }

  } # end for loop in cc



  ##########################################
  ## THEN get an estimated elementwise (K-cluster)
  ## run fewer iterations because takes longer to do Dahl calculations
  ##########################################

  All.Stuff.1 <- All.Stuff

  ###########

  All.Stuff <- All.Stuff.1
  All.Stuff.2 <- NULL

  ###########

  new.n.reps <- n.reps/2
  new.n.burn <- n.burn/2
  new.text <- paste("ROUND_2_",text,sep="")

  #
  All.Stuff.2$tau_0.v <- All.Stuff.2$tau.v <- All.Stuff.2$tau_int.v <- All.Stuff.2$G.v <- All.Stuff.2$K.v <- array(,new.n.reps)
  All.Stuff.2$row.flip.v  <- array(0,new.n.reps)

  All.Stuff.2$merge.flip.v  <- All.Stuff.2$split.changed.v <- array(0,new.n.reps)

  parm <- fn.init(true, data, zero.hemming, entropy.prop, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, col.DP.flag=FALSE, dahl=All.Stuff$dahl, d=mean(All.Stuff$d.v), max.d=max.d)
  init.parm <- parm

  err <- fn.quality.check(parm)
  if (err > 0)
    {stop(paste("failed QC at fn.init: err=",err))
    }


  flip.count <- 0


  All.Stuff.2$dahl <- NULL
  #### K.matrix <- array(0,c(init.parm$N,init.parm$N))
  r.matrix <- array(0,c(new.n.reps,init.parm$N))

  if (new.n.burn > 0)
  {
    for (cc in 1:new.n.burn)
    {parm <- fn.iter(data, parm, zero.hemming, entropy.prop, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, col.DP.flag=FALSE, split.merge.flag, max.d=max.d)


    if (cc %% 10 == 0)
      {print(paste(new.text, "BURN = ",cc,date(),"***********"))
      }
      }
    }


  for (cc in 1:new.n.reps)
    {parm <- fn.iter(data, parm, zero.hemming, entropy.prop, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, col.DP.flag=FALSE, split.merge.flag, max.d)

    All.Stuff.2$G.v[cc] <- parm$clust$G
    All.Stuff.2$K.v[cc] <- parm$clust$K
    All.Stuff.2$tau.v[cc] <- parm$tau
    All.Stuff.2$tau_0.v[cc] <- parm$tau_0
    All.Stuff.2$tau_int.v[cc] <- parm$tau_int

    # summarizing elementwise DP in "fn.groupwise.updates"

    All.Stuff.2$row.flip.v[cc]  <- parm$clust$row.flip

    r.matrix[cc,] <- parm$clust$s.v

    ### exact Dahl calculations: too expensive
    if (FALSE)
    {
    tmp.mat <- array(0,c(parm$N,parm$N))

    for (jj in 0:parm$clust$K)
      {indx.jj <- which(parm$clust$s.v==jj)
      tmp.mat[indx.jj,indx.jj] <- 1
      }

    K.matrix <- K.matrix + tmp.mat
    }

    if (cc %% 10 == 0)
      {print(paste(new.text, "REPS = ",cc,date(),"***********"))
    }

    } # end for loop in cc


  All.Stuff.2$parm <- parm
  All.Stuff.2$init.parm <- init.parm

  # normalize Dahl matrix
  #### K.matrix <- K.matrix/new.n.reps

  All.Stuff.2$dahl$min <- Inf

  # sample only up to 20 replications for Dahl calculations
  for (cc in sort(sample(new.n.reps,size=min(new.n.reps,20))))
  {

    ### exact Dahl calculations: too expensive
    if (FALSE)
    {
      tmp.mat <-array(0,c(parm$N,parm$N))

      for (jj in 0:All.Stuff.2$K.v[cc])
      {indx.jj <- which(r.matrix[cc,]==jj)
      tmp.mat[indx.jj,indx.jj] <- 1
      }

      dahl.dist <- sum((K.matrix - tmp.mat)^2)
      if (dahl.dist < All.Stuff.2$dahl$min)
      {All.Stuff.2$dahl$min <- dahl.dist
      All.Stuff.2$dahl$K.clust <- r.matrix[cc,]
      All.Stuff.2$dahl$K <- All.Stuff.2$K.v[cc]
      }
    }

    ### approx Dahl calculations
    if (TRUE)
    {
      dahl.dist <- 0
      abort.flag <- FALSE

      for (k in 1:max(r.matrix[cc,]))
      {i.k <- sort(which(r.matrix[cc,]==k))

      if (length(i.k) > 0)
      {for (i in i.k)
        for (j in i.k)
        {if (j > i)
        {K.ij <- mean(r.matrix[,i]==r.matrix[,j])
        delta.ij <- as.numeric(r.matrix[cc,i] == r.matrix[cc,i])
        dahl.dist <- dahl.dist + 2*(K.ij-delta.ij)^2
        if (dahl.dist > All.Stuff.2$dahl$min)
        {abort.flag <- TRUE
        break
        }
        }
        }
        if (abort.flag) break
      }
      }

      if (dahl.dist < All.Stuff.2$dahl$min)
      {All.Stuff.2$dahl$min <- dahl.dist
      All.Stuff.2$dahl$K.clust <- r.matrix[cc,]
      All.Stuff.2$dahl$K <- All.Stuff.2$K.v[cc]
      }
    }


  } # end for loop in cc


  # just for the record, although All.Stuff.1 and All.Stuff.2
  # have $d fixed at mean
  All.Stuff.1$d.v <- All.Stuff.2$d.v <- All.Stuff$d.v

  list(All.Stuff.1, All.Stuff.2)
  }
