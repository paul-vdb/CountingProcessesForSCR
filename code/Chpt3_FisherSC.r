########################
# Chapter 3 Examples
# Fisher Camera Trap Data
# Spatial Count Model with Data Augmentation, Reversible Jump, and Birth Death
########################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(ggplot2)
library(sf)
library(secr)
library(here)

setwd(here())

# Load all Nimble functions:
files <- dir("functions")
files <- grep('\\.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0("functions/", x))})
source(paste0("functions", "/", files[6]))

load("data/Fisher/FisherData.Rda")

buffer <- 6
traps <- fisher.data[["traps"]]

avg.dists <- apply(traps, 1, FUN = function(x){sqrt((traps[,1] - x[1])^2 + (traps[,2]-x[2])^2)})
avg.min <- apply(avg.dists, 1, FUN = function(x){min(x[x!=0])})
avg <- mean(avg.min)

# traps <- traps[-(1:13),]
xlim <- range(traps[,1])+c(-buffer,buffer)
ylim <- range(traps[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100

# Start: January 20, 2016
start.date <- as.Date("01-01-2016", format = "%d-%m-%Y")
end.date <- as.Date("05-03-2016", format = "%d-%m-%Y")

# 207 detections.
obs <- fisher.data[["observations"]]
obs <- obs[as.Date(obs$DateTime) >= start.date & end.date > as.Date(obs$DateTime),]
nrow(obs)

omega <- obs$TrapNumber
StudyPeriod <- as.numeric(end.date - start.date)

M <- 400
J <- nrow(traps)

# Chandler and Royle Spatial Count model, Algorithm 2
# Repeat of Burgar et al. 2018 analysis for paper.
counts <- numeric(J)
for(j in 1:J)
{
	counts[j] <- sum(obs$TrapNumber == j)
}

traps.secr <- read.traps(data= traps)
mask <- make.mask(traps, buffer = buffer)
a <- attr(mask, 'area')
nmask <- nrow(mask)
dists2 <- t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2}))

# Initialize model:
#-------------------------
inits <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	psi <- rbeta(1,1,1)
	z <- rbinom(M, prob=psi, size = 1)
	z[1:nrow(traps)] <- 1
	X[1:nrow(traps),] <- as.matrix(traps) + cbind(rnorm(nrow(traps), 0, 0.1), rnorm(nrow(traps), 0, 0.1))
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z
    )
}

# Run the same model from Burgar et al. Spatial Count on fisher.
#----------------------------------------------------------------------
constants <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
	area = area
	)

data <- list(
	z =  rep(NA, M),
	counts = counts
	)

## Data Augmentation SC Model:
SC_DA <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 50)
	lambda ~ dunif(0, 20)
	psi ~ dbeta(1,1)
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau2 )*z[k]
	}

	### Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time
		counts[j] ~ dpois(Hj[j])
	}
	
	## Just tracking the log-likelihood for my own reasons.
	ll <- sum(-Hj[1:J]) + sum(counts[1:J]*log(Hj[1:J]))
	N <- sum(z[1:M])
	D <- N/area
} )

Rmodel <- nimbleModel(SC_DA, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE) #, control = list(adaptive = FALSE, scaleWidth = 0.25))	

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M)	{
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE)	# Adaptive is okay for this model.
	# conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		# type = 'sampler_mys_grid', control = list(J = 64, nmask = nrow(mask), 
		# Time = StudyPeriod, mask = as.matrix(mask), count = counts, 
		# dist_squared = dists2))		
}

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(10000)
# mvSamples.da <- Cmcmc$mvSamples
# samples.da <- as.matrix(mvSamples.da)
# out.da <- mcmc(samples.da) #[-(1:5000),])
# plot(out.da[,c('N', 'D')])
# dev.new()
# max(out.da[,'ll'])
# summary(out.da)
# save(out.da, file = "geneticFisherMCMC15KXZ.Rda")

t.da <- system.time(
		outDA <- runMCMC(Cmcmc, niter = 100000, nburnin = 20000, nchains = 3, 
		thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)
# user  	system  elapsed 
# 8260.61   0.39    8261.13

# save(outDA, file = "../output/cameraFisherMCMC_DA.Rda")
effectiveSize(outDA)
efficiency.da <- effectiveSize(outDA)/t.da[1]
efficiency.da
# D         N          lambda     psi        sigma 
#0.05119420 0.05119420 0.31967348 0.05417244 0.08799885


# samplesDA <- runMCMC(Cmcmc, niter = 100000, nburnin = 20000, nchains = 3, 
	# thin = 1, inits = list(inits(), inits(), inits()))
# outDA <- mcmc.list(list(as.mcmc(samplesDA[[1]]), as.mcmc(samplesDA[[2]]), as.mcmc(samplesDA[[3]])))
# plot(outDA)
# save(outDA, file = "cameraFisherMCMC_DA.Rda")

################################
# Part 2) RJMCMC
################################

SC_RJ <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 50)
	lambda ~ dunif(0, 20)
	psi <- 0.5	## Trick to remove prior on N while still using indicator variable...
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau2 )*z[k]
	}

	### Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time
		counts[j] ~ dpois(Hj[j])
	}
	ll <- sum(-Hj[1:J]) + sum(counts[1:J]*log(Hj[1:J]))
	N <- sum(z[1:M])
	D <- N/area
} )

Rmodel <- nimbleModel(SC_RJ, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'N', 'D'))

conf$removeSamplers('X')
for(i in 1:M)	{
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'sampler_myzs_RJMCMC', control = list(scale = 1.7))
}
conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE) #, control = list(adaptive = FALSE, scaleWidth = 0.25))	

## Need to pass the prior limits, name of activity centre used and then the +/- tuning parameter delta to 
## to the RJMCMC sampler.
conf$removeSamplers('z')
conf$addSampler(target = 'z[1:400]', type = 'sampler_myRJMCMC_z', 
		control = list(ActivityCentre = "X", delta = 7, xlim = xlim, ylim = ylim))	

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(30000, time = TRUE)
# times.rjmcmc <- Cmcmc$getTimes()
# mvSamples <- Cmcmc$mvSamples
# samples.rj <- as.matrix(mvSamples)[-(1:5000),]
# out.rj <- mcmc(samples.rj)
# plot(out.rj[,c('sigma','lambda')])
# plot(out.rj[,c('N', 'D')])

t.rj <- system.time(
		outRJ <- runMCMC(Cmcmc, niter = 100000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)
save(outRJ, file = "../output/cameraFisherMCMC_RJ.Rda")

# user  	system  elapsed 
# 3626.34   0.71    3632.03 

effectiveSize(outRJ)
efficiency.rj <- effectiveSize(outRJ)/t.rj[1]
efficiency.rj
plot(outRJ[,'N'])


################################
# Part 3) Birth Death
################################

SC_BD <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 50)
	lambda ~ dunif(0, 20)
	psi <- 0.5	## Trick to remove prior on N while still using indicator variable...
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		Hkj[k,1:J] <- lambda*exp(-d2[k,1:J] *tau2 )*Time
		Hkjz[k,1:J] <- Hkj[k,1:J]*z[k]
	}

	### Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(Hkjz[1:M,j])
		counts[j] ~ dpois(Hj[j])
	}
	N <- sum(z[1:M])
	D <- N/area
} )


Rmodel <- nimbleModel(SC_BD, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'N', 'D'))

## Custom sampler than immediately returns if z = 0 saving a lot of time.
## Could speed up with a joint sampler but I think c++ loops are really fast if they are empty so
## why bother.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'sampler_myzs_BD', control = list(scale = 1.7, xlim=xlim, ylim=ylim))
}
conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE) #, control = list(adaptive = FALSE, scaleWidth = 0.25))	

## Need to pass the prior limits, name of activity centre used and then the +/- tuning parameter delta to 
## to the RJMCMC sampler.
conf$removeSamplers('z')
conf$addSampler(target = 'z[1:400]', type = 'sampler_myBDSC_z',
		control = list(t0 = 1, BirthRate = 2.5, M = M, ObservedCounts = counts))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(50000, time = TRUE)
# Cmcmc$getTimes()
# mvSamples.bd <- Cmcmc$mvSamples
# samples.bd <- as.matrix(mvSamples.bd)[-(1:5000),]
# out.bd <- mcmc(samples.bd)
# plot(out.bd[,c('sigma','lambda')])
# plot(out.bd[,c('N', 'D')])
# plot(out.bd[,'ll'])
# summary(out.bd[,"N"])

# samplesBD <- runMCMC(Cmcmc, niter = 100000, nburnin = 20000, nchains = 3, 
	# thin = 1, inits = list(inits(), inits(), inits()))
# outBD <- mcmc.list(list(as.mcmc(samplesBD[[1]]), as.mcmc(samplesBD[[2]]), as.mcmc(samplesBD[[3]])))
# plot(outBD)

t.bd <- system.time(
		outBD <- runMCMC(Cmcmc, niter = 100000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
)
# save(outBD, file = "../output/cameraFisherMCMC_BD.Rda")

# user  	system  elapsed 
# 7069.47   0.48    7070.16

effectiveSize(outBD)
efficiency.bd <- effectiveSize(outBD)/t.bd[1]
efficiency.bd
plot(outBD[, "N"])