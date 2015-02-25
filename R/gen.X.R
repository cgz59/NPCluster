
gen.X <- function(n, p, true_parm) {
	
	data <- NULL
	data$X <- array(,c(n,p))
	
	for (xx in 1:p)
	{s.v <- true_parm$clust$s.mt[,true_parm$clust$c.v[xx]]
	 
	 x.v <- array(,n)
	 indx.0 <- which(s.v==0)
	 indx.pos <- which(s.v>0)
	 
	 if (length(indx.0)>0)
	 {x.v[indx.0] <- rnorm(n=length(indx.0),sd=true_parm$tau_0)
	 }
	 
	 x.v[indx.pos] <- rnorm(n=length(indx.pos),mean=true_parm$clust$phi.v[s.v[indx.pos]],sd=true_parm$tau)
	 
	 data$X[,xx] <- x.v	
	}
	
	
	
	###########################################
	# small proportion of missing X values
	###########################################
	
	data$X.missing.x <- NULL
	data$X.missing.y <- NULL
	
	###########################################
	# random split of 100 X prop % missing
	###########################################
	
	n2 <- n
	n1 <- 0
	data$missing.indx <- NULL
	data$non.missing.indx <- 1:n
	num.X.miss <- 0
	
	K.max <- round(n2/2)
	G.max <- round(p/2) 
	
	###########################################
	# dummy responses
	###########################################
	
	data$Y <- rep(0,n2)
	data$delta <- rep(0,n2)
	data$true <- NULL
	data$true$Y <- data$Y
	data$true$delta <- data$delta
	
	############
	
	true <- NULL
	true$a.R <- true_parm$clust$M
	true$b0 <- 2.2
	true$b1 <- true_parm$b1
	
	#########################################
	# generating the R- and C- clusters
	########################################
	
	true$shift <- 1e-4
	
	return(list(data = data, true = true))
}
