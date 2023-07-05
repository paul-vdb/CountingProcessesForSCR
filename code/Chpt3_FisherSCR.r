########################
# Chapter 3 Examples
# Fisher Hair Snag Data
# Spatial Capture-Recapture Model with Data Augmentation, Reversible Jump, and Birth Death
########################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(ggplot2)
library(sf)
# devtools::install_github("nimble-dev/compareMCMCs", subdir = "compareMCMCs")
library(compareMCMCs)

library(here)
setwd(here())

# Load Nimble Samplers:
files <- dir("functions")
files <- grep('\\.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0("functions/", x))})
source(paste0("functions", "/", files[6]))

tdf <- read.csv("data/Fisher/tdf.csv", header = TRUE, sep = ",", check.names = FALSE) ## NEED CHECK.NAMES = FALSE otherwise 3d will not work
edf <- read.csv("data/Fisher/edf.csv", header = TRUE, sep = ",", check.names = FALSE) ## NEED CHECK.NAMES = FALSE otherwise 3d will not work
edf <- edf[order(edf$individualID_2016),]
animal <- edf[!duplicated(edf$individualID_2016),]
sex <- animal$sex
sex <- ifelse(sex=="M",0,1) # 0 is female, 1 is male
# Now put it on the animal:
sex[unique(edf$individualID_2016)]
traps <- tdf[,2:3]
coord.scale <- 1000 # 1 km units
buffer <- 6 #6 km unit buffer - same as manuscript...
traps <- traps/coord.scale

K <- max(edf$ind)

xlim = range(traps[,1])+c(-buffer,buffer)
ylim = range(traps[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100	# Density reported per 100 sq km

M <- 400
y <- array(0, dim = c(M, 64, 4))
# Add the captures as 1s with the good old cbind trick.
y[cbind(edf$individualID_2016,edf$secr_trapID, edf$occasionID_2016)] <- 1
sum(y) == nrow(edf) 	# It worked right?

y_sum <- apply(y, 1:2, sum)

################################
# Part 1) Data Augmentation
################################

inits <- function()
{
	psi <- rbeta(1, 1+K, M-K)
	z <- numeric(M)
	z[1:K] <- NA
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
	for( i in 1:K ){
		X[i,] <- as.numeric(traps[which.max(y_sum[i,]),] + rnorm(2, 0, 0.1))
	}
	list(X = X, z = z, sigma = runif(1,1,3), lambda = runif(1,0.05,0.2), psi = psi)
}


SCR_DA <- nimbleCode({
	sigma ~ dunif(0,1000) # uninformative prior
	psi ~ dbeta(1,1)
	lambda ~ dunif(0,20)
	tau2 <- 1/(2*sigma^2)
	
	for(i in 1:M){
		z[i] ~ dbern(psi)
		X[i,1]~dunif(xlim[1],xlim[2])
		X[i,2]~dunif(ylim[1],ylim[2])
		d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		Hkj[i,1:J] <- 30*lambda*exp(-d2[i, 1:J]*tau2)
		pkj[i,1:J]<- z[i]*(1-exp(-Hkj[i,1:J]))
	
		#From Daniel Turek in nimbleSCR package. Fast binomial! Avoids loopin.
		y[i, 1:J] ~ dbinom_vector(size = trials[1:J], prob = pkj[i,1:J])
	}
	## Impute counts for comparison:
	N <- sum(z[1:M])
	D <- N/area
})

constants <- list(
    J = nrow(traps),
	area = area,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    M = M,
	trials = rep(4, nrow(traps))
	)

K <- max(edf$individualID_2016)
data <- list(
	z =  c(rep(1, K), rep(NA, M-K)),
	y = y_sum
	)

Rmodel <- nimbleModel(SCR_DA, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M){
	if(i <= K){	conf$addSampler(target = paste0('X[', i, ', 1:2]'),
			type = 'RW_block', silent = TRUE)
	}else{
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'RW_block', control = list(adapt = FALSE, scale = 1))
	}
}

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

t.da <- system.time(
	 outDA <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, 
		thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)

# user  	system  elapsed 
# 6763.38   2.72    6777.28

# save(outDA, file = "../output/geneticFisherMCMC_DA.Rda")
effectiveSize(outDA)
efficiency.d <- effectiveSize(outDA)/t.da[1]
efficiency.d
# D         N         lambda    psi       sigma 
# 3.0743356 3.0743356 1.0137280 3.3875988 0.7036175
plot(outDA[, "N"])

# Cmcmc$run(10000)
# mvSamples.da <- Cmcmc$mvSamples
# samples.da <- as.matrix(mvSamples.da)
# out.da <- mcmc(samples.da[-(1:5000),])
# plot(out.da[,c('N', 'D')])
# dev.new()
# plot(out.da[,c('sigma','lambda')])
# summary(out.da)
# save(out, file = "geneticFisherMCMC.Rda")

################################
# Part 2) RJMCMC
################################

SCR_RJMCMC <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 10)
	lambda ~ dunif(0, 10)
	psi <- 0.5	## Trick to remove prior on N while still using indicator variable...
	tau2 <- 1/(2*sigma^2)
	
	# Observed data likelihood
	for(i in 1:M){
		z[i] ~ dbern(psi)
		X[i,1]~dunif(xlim[1],xlim[2])
		X[i,2]~dunif(ylim[1],ylim[2])
		d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		
		Hkj[i,1:J] <- 30*lambda*exp(-d2[i, 1:J]*tau2)
		pkj[i,1:J]<- z[i]*(1-exp(-Hkj[i,1:J]))
	
		#From Daniel Turek in nimbleSCR package. Fast binomial! Avoids loopin.
		y[i, 1:J] ~ dbinom_vector(size = trials[1:J], prob = pkj[i,1:J])
	}

	N <- sum(z[1:M])
	D <- N/area
} )

constants <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    M = M,
	area = diff(xlim)*diff(ylim)/100,
	trials = rep(4, nrow(traps))
	)

data <- list(
	z =  c(rep(1, K), rep(NA, M-K)),
	y = y_sum
	)

inits <- function()
{
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	traps.det <- apply(y_sum[1:K,], 1, which.max)
	X[1:K,] <- as.matrix(traps[traps.det,]) + cbind(rnorm(K, 0, 0.5),rnorm(K, 0, 0.5))			  
	z <- rep(0, M)
	z[1:K] <- NA
	z[(K+1):(K+10)] <- 1
	list(
		lambda = lambda,
		sigma = sigma,
		X = X,
		z = z
    )	
}

Rmodel <- nimbleModel(SCR_RJMCMC , constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'N', 'D'))

## Custom sampler than immediately returns if z = 0 saving a lot of time.
## Could speed up with a joint sampler but I think c++ loops are really fast if they are empty so
## why bother.
conf$removeSamplers('X')
for(i in 1:M){
	if(i <= K){	conf$addSampler(target = paste0('X[', i, ', 1:2]'),
			type = 'RW_block', silent = TRUE)
	}else{
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'sampler_myzs_RJMCMC', control = list(scale = 1))
	}
}

## Need to pass the prior limits, name of activity centre used and then the +/- tuning parameter delta to 
## to the RJMCMC sampler.
conf$removeSamplers('z')
conf$addSampler(target = paste0('z[',(K+1),":",M,"]"), type = 'sampler_myRJMCMC_z',
		control = list(ActivityCentre = "X", delta = 5, xlim = xlim, ylim = ylim))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


t.rj <- system.time(
	 outRJ <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, 
		thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)

# user  	system  elapsed 
# 5106.58   0.95    5108.17

# save(outRJ, file = "../output/geneticFisherMCMC_RJ.Rda")
effectiveSize(outRJ)
efficiency.rj <- effectiveSize(outRJ)/t.rj[1]
efficiency.rj
# D         N         lambda    sigma 
# 0.1177652 0.1177652 0.8002378 0.8109321
plot(outRJ[, "N"])


# samplesRJMCMC <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, 
	# thin = 1, inits = list(inits(), inits(), inits()))
# outRJ <- mcmc.list(list(as.mcmc(samplesRJMCMC[[1]]), as.mcmc(samplesRJMCMC[[2]]), as.mcmc(samplesRJMCMC[[3]])))
# plot(outRJ)

# Cmcmc$run(10000, time = TRUE)
# Cmcmc$getTimes()
# mvSamples.rjmcmc <- Cmcmc$mvSamples
# samples.rjmcmc <- as.matrix(mvSamples.rjmcmc)[-(1:5000),]
# out.rjmcmc <- mcmc(samples.rjmcmc)
# plot(out.rjmcmc[,c('sigma','lambda')])
# plot(out.rjmcmc[,c('N', 'D')])


################################
# Part 3) Birth Death
################################

SCR_BD <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 10)
	lambda ~ dunif(0, 10)
	psi <- 0.5	## Trick to remove prior on N while still using indicator variable...
	tau2 <- 1/(2*sigma^2)
	
	# Observed data likelihood
	for(i in 1:M){
		z[i] ~ dbern(psi)
		X[i,1]~dunif(xlim[1],xlim[2])
		X[i,2]~dunif(ylim[1],ylim[2])
		d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		
		Hkj[i,1:J] <- 30*lambda*exp(-d2[i, 1:J]*tau2)
		Hk[i] <- sum(Hkj[i,1:J])*4				# Needed this hard coded (w/ this name) for the internal z sampler, 4 repeated trials!!
		pkj[i,1:J]<- z[i]*(1-exp(-Hkj[i,1:J]))
	
		#From Daniel Turek in nimbleSCR package. Fast binomial! Avoids loopin.
		y[i, 1:J] ~ dbinom_vector(size = trials[1:J], prob = pkj[i,1:J])
	}
	
	N <- sum(z[1:M])
	D <- N/area
} )


constants <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    M = M,
	area = diff(xlim)*diff(ylim)/100,
	trials = rep(4, nrow(traps))
	)

data <- list(
	z =  c(rep(1, K), rep(NA, M-K)),
	y = y_sum
	)

inits <- function()
{
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	traps.det <- apply(y_sum[1:K,], 1, which.max)
	X[1:K,] <- as.matrix(traps[traps.det,]) + cbind(rnorm(K, 0, 0.5),rnorm(K, 0, 0.5))
	z <- rep(0, M)
	z[1:K] <- NA
	z[(K+1):(K+10)] <- 1
	list(
		lambda = lambda,
		sigma = sigma,
		X = X,
		z = z
    )	
}

Rmodel <- nimbleModel(SCR_BD, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

## Custom sampler than immediately returns if z = 0 saving a lot of time.
## Could speed up with a joint sampler but I think c++ loops are really fast if they are empty so
## why bother.
conf$removeSamplers('X')
for(i in 1:M){
	if(i <= K){	conf$addSampler(target = paste0('X[', i, ', 1:2]'),
			type = 'RW_block', silent = TRUE)
	}else{
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'sampler_myzs_BD', control = list(scale = 1, xlim=xlim, ylim=ylim))
	}
}

## Need to pass the prior limits, name of activity centre used and then the +/- tuning parameter delta to 
## to the RJMCMC sampler.
conf$removeSamplers('z')
conf$addSampler(target = paste0('z[',(K+1),":",M,"]"), type = 'sampler_myBD_z',
		control = list(t0 = 1, BirthRate = 2, ObservedZ = c(rep(1,K), rep(0, M-K))))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

t.bd <- system.time(
	 outBD <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, 
		thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)

# user  	system  elapsed 
# 7290.06   2.70    7297.11

# save(outBD, file = "../output/geneticFisherMCMC_BD.Rda")
effectiveSize(outBD)
efficiency.bd <- effectiveSize(outBD)/t.bd[1]
efficiency.bd
# D         N         lambda    sigma 
# 0.1177652 0.1177652 0.8002378 0.8109321
plot(outBD[, "N"])
summary(outBD)

samplesBD <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()))
outBD <- mcmc.list(list(as.mcmc(samplesBD[[1]]), as.mcmc(samplesBD[[2]]), as.mcmc(samplesBD[[3]])))
# plot(outBD)
# save(outBD, file = "../output/geneticFisherMCMC_BD.Rda")

# load("geneticFisherMCMC_BD3.Rda")
# outBD3 <- outBD
# load("geneticFisherMCMC_BD.Rda")

# Cmcmc$run(5000, time = TRUE)
# Cmcmc$getTimes()
# mvSamples.bd <- Cmcmc$mvSamples
# samples.bd <- as.matrix(mvSamples.bd) #[-(1:5000),]
# out.bd <- mcmc(samples.bd)
# plot(out.bd[,c('sigma','lambda')])
# plot(out.bd[,c('N', 'D')])
