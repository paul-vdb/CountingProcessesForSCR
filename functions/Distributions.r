##--------------------------------------------------------------------------
## Dummy distribution for activity centres.
##--------------------------------------------------------------------------
dX <- nimbleFunction(
    run = function(x = double(1), log = integer(0, default = 0)) {
        returnType(double(0))
        return(0)
    }
)
##--------------------------------------------------------------------------
## Random Generator not implemented
##--------------------------------------------------------------------------
rX <- nimbleFunction(
    run = function(n = integer(0)) {
        print('Error: rX() function not implemented')
        returnType(double(1))
        return(c(0, 0))
    }
)

# Stupid holder function for sampling ID.
# I've added pID here so that I can put the z dependency here as a prior on the
# ID and make the categorical sampler a tiny bit quicker. Can also
# include other hard prior info here making it helpful.

##--------------------------------------------------------------------------
## Dummy distribution to allow different priors on a categorical variable without
## declaring it as such so we can cheat if it sums to 1.
##--------------------------------------------------------------------------
dID <- nimbleFunction(
    run = function(x = integer(0), pID = double(1), log = integer(0, default = 0)) {
        returnType(double(0))
		if(log) return(log(pID[x])) else return(pID[x])
    }
)
##--------------------------------------------------------------------------
## Simple categorical random generator.
##--------------------------------------------------------------------------
rID <- nimbleFunction(
    run = function(n = integer(0), pID = double(1)) {
        returnType(integer(0))
        return(rcat(1, pID/max(pID)))
    }
)

##--------------------------------------------------------------------------
# Dummy distribution but avoids the problem of using the ones trick 
# and having to divide by a constant. Not sure impact on speed
# but all it does is index the detection probability for a given trap number.
##--------------------------------------------------------------------------
dTrap <- nimbleFunction(
    run = function(x = integer(0), p = double(2), ID = integer(0), log = integer(0, default = 0)) {
        returnType(double(0))
		if(log) return(log(p[ID, x])) else return(p[ID, x])
    }
)


##--------------------------------------------------------------------------
## Didn't really implement this.
##--------------------------------------------------------------------------
rTrap <- nimbleFunction(
    run = function(n = integer(0), p = double(2), ID = integer(0)) {
        returnType(integer(0))
        return(rcat(1, p[ID,]/max(p[ID,])))
    }
)


##--------------------------------------------------------------------------
## Dirty little internal trick to add whatever probs anyway I please and avoiding the zero or ones trick.
##--------------------------------------------------------------------------
dTrick <- nimbleFunction(
    run = function(x = integer(0), p = double(0), log = integer(0, default = 0)) {
        returnType(double(0))
		if(log) return(log(p)) else return(p)
    }
)

##--------------------------------------------------------------------------
## Return something....
##--------------------------------------------------------------------------
rTrick <- nimbleFunction(
    run = function(n = integer(0), p = double(0)) {
        returnType(integer(0))
        return(rexp(1, p))
    }
)



##--------------------------------------------------------------------------
## Distribution with cue production time marginalized:
##--------------------------------------------------------------------------
dnorm_vector_marg <- nimbleFunction(
  run = function( x = double(1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
	m <- sum(y)
	if(m == 1){
		logProb <- 0
	}else{
		td <- (x - mean)*y
		etd <- sum(td)/m
		logProb <- (1-m)*log(sd) - sum(y*(td - etd)^2)/(2*sd^2)
	}	
    if(log) return(logProb) else return(exp(logProb))
  })

##--------------------------------------------------------------------------
## Randomization with cue production time marginalized.
## Not fully implemented
##--------------------------------------------------------------------------
rnorm_vector_marg <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1)
  ) {
    returnType(double(1))
    return(rnorm(length(y), mean = mean, sd = sd))
  })

##--------------------------------------------------------------------------
## Distribution to vectorize the capture history conditional on being detected
##--------------------------------------------------------------------------
dbinom_vector_ascr <- nimbleFunction(
  run = function( x = double(1),
                  size = double(1),
                  prob = double(1), 
				  pcapt = double(0),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE)) - log(pcapt)
    if(log) return(logProb) else return(exp(logProb))
  })

##--------------------------------------------------------------------------
## Random distribution for conditional bernoulli
##--------------------------------------------------------------------------
rbinom_vector_ascr <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  size = double(1),
                  prob = double(1),
				  pcapt = double(0)
  ) {
    returnType(double(1))
	capthist <- rbinom(length(size), prob = prob, size = size)
	capthist[rcat(1, prob)] <- 1	# Ensures at least one detection.
    return(capthist)
  })


##--------------------------------------------------------------------------
## Adding the Poisson-Binomial distribution to Nimble from R.
##--------------------------------------------------------------------------
RdensityFunction <- function(x = double(1), prob = double(2), z = double(1), log = integer(0, default = 0)) {
    require(poisbinom)
	val <- 0
	for(j in 1:ncol(prob))
	{
		val <- val + dpoisbinom(x[j], pp = prob[z==1,j], log = 1)
    }
	returnType = double(0)
    if(log) return(val) else return(exp(val))
}

dPoisBin <- nimbleRcall(
    prototype = function(x = double(1), prob = double(2), z = double(1), log = integer(0)) {},
    returnType = double(0),
    Rfun = 'RdensityFunction'
)

##--------------------------------------------------------------------------
## Distribution to speed up Spatial Count model. Didn't help much.
##--------------------------------------------------------------------------
dPoisSC <- nimbleFunction(
  run = function( x = double(1),
                  lambda = double(2),
				  J = integer(0),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
	val <- 0
	for(j in 1:J)
	{
		Lamj <- sum(lambda[,j])
		val <- val + x[j]*log(Lamj) - Lamj
    }
    if(log) return(val) else return(exp(val))
  })

##--------------------------------------------------------------------------
## Vectorized random Poisson.
##--------------------------------------------------------------------------
rPoisSC <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  lambda = double(2), J = integer(0)) {
    returnType(double(1))
	nj <- numeric(J)
	for(j in 1:J)
	{
		nj[j] <- rpois(1, sum(lambda[,j]))
	}
    return(nj)
  })


##--------------------------------------------------------------------------
## Distribution to vectorize the Poisson counts
##--------------------------------------------------------------------------
dpois_vector <- nimbleFunction(
  run = function( x = double(1),
                  lambda = double(1),
				  J = double(0),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
    logProb <- sum(dpois(x, lambda, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })

##--------------------------------------------------------------------------
## Random distribution for vector Poisson
##--------------------------------------------------------------------------
rpois_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  lambda = double(1),
				  J = double(0))
  {
    returnType(double(1))
    return(rpois(J, lambda))
  })



##--------------------------------------------------------------------------
## Register Distributions
##--------------------------------------------------------------------------
registerDistributions(
    list(dX = list(BUGSdist = 'dX()',
                   types = c('value = double(1)')),
		 dID = list(BUGSdist = 'dID(pID)',
                   types = c('value = integer(0)', 'pID = double(1)')),
		 dTrap = list(BUGSdist = 'dTrap(p, ID)',
                   types = c('value = integer(0)', 'p = double(2)', 'ID = integer(0)')),
		 dTrick = list(BUGSdist = 'dTrick(p)',
                   types = c('value = integer(0)', 'p = double(0)')),			   
		 dnorm_vector_marg = list(BUGSdist = 'dnorm_vector_marg(mean, sd, y)',
				   types = c('value = double(1)', 'mean = double(1)', 'sd = double(0)', 'y = double(1)')),
      	 dPoisBin = list(
				   BUGSdist = 'dPoisBin(prob, z)', Rdist = 'dPoisBin(prob, z)',
				   types = c('value = double(1)', 'prob = double(2)', 'z = double(1)')
				   ),
      	 dPoisSC = list(
				   BUGSdist = 'dPoisSC(lambda,J)', Rdist = 'dPoisSC(lambda,J)',
				   types = c('value = double(1)', 'lambda = double(2)', 'J = integer(0)')
				   ),
		dbinom_vector_ascr = list(
				   BUGSdist = 'dbinom_vector_ascr(size, prob, pcapt)', Rdist = 'dbinom_vector_ascr(size, prob, pcapt)',
				   types = c('value = double(1)','size = double(1)','prob = double(1)', 'pcapt = double(0)')
				   ),
		dpois_vector = list(
				   BUGSdist = 'dpois_vector(lambda, J)', Rdist = 'dpois_vector(lambda, J)',
				   types = c('value = double(1)', 'lambda = double(1)', 'J = double(0)')
				   )		   
		)
)