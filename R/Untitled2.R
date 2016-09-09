PDP_fn.log.lik <- function(gg, x.mt, parm, colSums = TRUE)
{
  if (gg > 0)
  {
    a2.v <- parm$clust$A.mt[,gg]
    z.g.v <- parm$clust$s.mt[,gg] > 0
    log.lik.v <- rep(0,ncol(x.mt))

    if (sum(z.g.v) > 0)
    {
      a2.1.v <- parm$clust$A.mt[z.g.v,gg]
      small.X.1 <- matrix(x.mt[z.g.v,], ncol = ncol(x.mt)) # HOT

      if (colSums)
      {
        if (!parm$flip.sign)
        {
          tmp <- colSums(-.5 * (small.X.1 - a2.1.v) ^ 2)
        }
        if (parm$flip.sign)
        {
          tmp <-
            colSums(-.5 * (small.X.1 - a2.1.v) ^ 2) + colSums(-.5 * (-small.X.1 - a2.1.v) ^
                                                                2)
        }
      } # end 	if (colSums)

      if (!colSums)
      {
        if (!parm$flip.sign)
        {
          tmp <- sum(-.5 * (small.X.1 - a2.1.v) ^ 2)
        }
        if (parm$flip.sign)
        {
          tmp <-
            sum(-.5 * (small.X.1 - a2.1.v) ^ 2) + sum(-.5 * (-small.X.1 - a2.1.v) ^
                                                        2)
        }
      } # end if (!colSums)

      log.lik.v <-
        log.lik.v + tmp / (parm$tau ^ 2) + sum(z.g.v) * (-.5 * log(2 * pi) - log(parm$tau))
    } # end if (sum(z.g.v) > 0)

    if (sum(1 - z.g.v) > 0)
    {
      small.X.0 <- matrix(x.mt[!z.g.v,], ncol = ncol(x.mt))
      if (colSums)
      {
        if (!parm$flip.sign)
        {
          tmp <- colSums(-.5 * small.X.0 ^ 2)
        }
        if (parm$flip.sign)
        {
          tmp <- 2 * colSums(-.5 * small.X.0 ^ 2)
        }
      }
      if (!colSums)
      {
        if (!parm$flip.sign)
        {
          tmp <- sum(-.5 * small.X.0 ^ 2)
        }
        if (parm$flip.sign)
        {
          tmp <- 2 * sum(-.5 * small.X.0 ^ 2)
        }
      }
      log.lik.v <-
        log.lik.v + tmp / (parm$tau_0 ^ 2) + sum(1 - z.g.v) * (-.5 * log(2 * pi) -
                                                               log(parm$tau_0))
    } # end if (sum(1-z.g.v) > 0)
  } # end if (gg > 0)

  if (gg == 0)
  {
    a2.v <- rep(1,parm$n2)
    small.X <- x.mt
    if (colSums)
    {
      if (!parm$flip.sign)
      {
        tmp <- colSums(-.5 * (small.X - a2.v) ^ 2)
        save <- tmp
      }
      if (parm$flip.sign)
      {
        tmp <-
          colSums(-.5 * (small.X - a2.v) ^ 2) + colSums(-.5 * (-small.X - a2.v) ^
                                                          2)
      }
    }
    if (!colSums)
    {
      if (!parm$flip.sign)
      {
        tmp <- sum(-.5 * (small.X - a2.v) ^ 2)
      }
      if (parm$flip.sign)
      {
        tmp <- sum(-.5 * (small.X - a2.v) ^ 2) + sum(-.5 * (-small.X - a2.v) ^ 2)
      }
    } # end if (gg == 0)

    log.lik.v <-
      tmp / parm$tau_int ^ 2 - parm$n2 * .5 * log(2 * pi) - parm$n2 * log(parm$tau_int)
  }

  log.lik.v
}
