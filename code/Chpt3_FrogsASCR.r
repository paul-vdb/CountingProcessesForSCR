########################
# Chapter 3 Examples
# Frog Acoustic Data
# Acoustic Spatial Capture-Recapture Model comparing DA, Conditional, and Semi-Complete
########################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)
library(ggplot2)
library(secr)
library(ascr)
library(MCMCvis)
# devtools::install_github("nimble-dev/compareMCMCs", subdir = "compareMCMCs")
library(compareMCMCs)

library(here)
setwd(here())

# Load Nimble Samplers:
files <- dir("functions")
files <- grep('\\.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0("functions/", x))})
source(paste0("functions", "/", files[6]))

load("data/Frog/stacked-lightfooti.Rdata")
# traps2 <- convert.traps(traps)
# mask <- make.mask(traps = traps2, buffer = 15, spacing = 1, type = "trapbuffer")


#########################################
## Nimble Model Definitions:
#########################################

## Conditional ID-ASCR model
ascr_id_cond <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
	g0 ~ dunif(0, 50)
	
    for(k in 1:K) {
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])

		## Distances
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		expTime[k,1:J] <- sqrt(d2[k,1:J])/nu

		## Prob detected
        pkj[k,1:J] <- (1-exp(-g0*exp(-d2[k,1:J]*tau2)))

        ## Prob detected at least once.
        pk.[k] <-(1-prod(1-pkj[k,1:J]))
		
		## Count process
		counts[k] ~ dpois(pk.[k]*lambda*Time)
    }
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector_ascr(size = trials[1:J], pkj[ID[i],1:J],  pcapt = pk.[ID[i]])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
    }

	# ESA <- ESA_A(d2mask[1:nmask, 1:J], lambda*Time, sigma, g0, area = A, detfn = 2, nmask)
	ESA <- ESA_A_CPP(d2mask[1:nmask, 1:J], lambda*Time, sigma, g0, area = A, detfn = 2)
	panimal <- ESA/area
	p0 <- ESA^(-K)
	zero ~ dTrick(p0)

    # Predicted population size
	# sdK <- sqrt(K*(1-panimal))
	# rK ~ dnorm(K, sd = sdK)
	# D <- rK/ESA
})

## Semi-Complete ID-ASCR model
ascr_id_SC <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
	g0 ~ dunif(0, 50)
	psi ~ dbeta(1,1)
	
    for(k in 1:K) {
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])

		## Distances
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		expTime[k,1:J] <- sqrt(d2[k,1:J])/nu

		## Prob detected
        pkj[k,1:J] <- (1-exp(-g0*exp(-d2[k,1:J]*tau2)))

        ## Hazard rate for animal across all traps.
        pk.[k] <-(1-prod(1-pkj[k,1:J]))
		
		## Count process
		counts[k] ~ dpois(pk.[k]*lambda*Time)
		
    }	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector_ascr(size = trials[1:J], pkj[ID[i],1:J],  pcapt = pk.[ID[i]])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
    }

	# ESA <- ESA_A(d2mask[1:nmask, 1:J], lambda*Time, sigma, g0, area = A, detfn = 2, nmask)
	ESA <- ESA_A_CPP(d2mask[1:nmask, 1:J], lambda*Time, sigma, g0, area = A, detfn = 2)
	panimal <- ESA/area
	p0 <- panimal^(-K)
	zero ~ dTrick(p0)

	## Unobserved Animals:
	# N0 ~ dnegbin(panimal, K)
	K ~ dbinom(size = 2*300, prob = psi*panimal)
	
    # Predicted population size
	D <- (psi*300)/area
})

## Data Augmented ID-ASCR model
ascr_id_DA <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    sigma ~ dunif(0, 10)  # Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
	g0 ~ dunif(0, 50)
	psi ~ dbeta(1,1)

    for(k in 1:M) {
		z[k] ~ dbern(psi)
	
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])

		## Distances
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		expTime[k,1:J] <- sqrt(d2[k,1:J])/nu

		## Prob detected
        pkj[k,1:J] <- (1-exp(-g0*exp(-d2[k,1:J]*tau2)))

        ## Hazard rate for animal across all traps.
        pk.[k] <-(1-prod(1-pkj[k,1:J]))*z[k]
		
		## Count process
		counts[k] ~ dpois(pk.[k]*lambda*Time)
		
    }	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector_ascr(size = trials[1:J], pkj[ID[i],1:J],  pcapt = pk.[ID[i]])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
    }

	N <- sum(z[1:M])
	
    # Predicted population size
	D <- (psi*M/2)/area
})

#########################################
## Data Input:
#########################################
## Let's make a square mask for our purposes of comparisons.
traps.secr <- read.traps(data = data.frame(traps))
mask <- make.mask(traps.secr, buffer = 15)
a <- attr(mask, 'area')
nmask <- nrow(mask)
d2mask <- t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2}))

delta <- sqrt(attr(mask, 'area'))/2 ## Correction factor for the bounds of centroid of mask points and actual limits.
xlim <- range(mask[,1]) + c(-delta, delta)
ylim <- range(mask[,2]) + c(-delta, delta)
area <- nrow(mask)*attr(mask, 'area')
toa <- capt.all$toa
capt <- capt.all$bincapt[,1:6]
# tmin <- apply(toa, 1, max)
# keep <- which(tmin < 1200)
# toa <- toa[keep,]
ID <- capt.all$bincapt[, 7]
# ID <- ID[keep]
ID <- as.integer(as.factor(ID))
IDBen <- ID
# capt <- capt[keep,]
K <- max(ID)

counts <- as.numeric(table(ID))

# Constants:
M <- 300*2
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30
n_obs = n

constants_ID <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
	n_obs = n_obs,
	trials = rep(1, J),
	nu = nu,
	area = area,
	A = attr(mask, 'area'),
	d2mask = d2mask,
	nmask = nmask,
	ID = ID,
	K = K
	)

constants_IDSC <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
	n_obs = n_obs,
	trials = rep(1, J),
	nu = nu,
	area = area,
	A = attr(mask, 'area'),
	d2mask = d2mask,
	nmask = nmask,
	ID = ID
	)

constants_IDDA <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
	n_obs = n_obs,
	trials = rep(1, J),
	nu = nu,
	area = diff(xlim)*diff(ylim)/10000,
	A = attr(mask, 'area'),
	d2mask = d2mask,
	nmask = nmask,
	ID = ID,
	K = K
	)
##-------------------

## Density Data
data_ID <- list(
	zero = 0,
    y = capt,
	toa = toa,
	counts = counts
)

data_IDSC <- list(
	zero = 0,
    y = capt,
	toa = toa,
	counts = counts,
	N0 = NA,
	K = K
)

data_IDDA <- list(
    y = capt,
	toa = toa,
	counts = c(counts, rep(0, M-K)),
	z = c(rep(1, K), rep(NA, M-K))
)
##-------------------

##########################################
## MCMC Initializations
##########################################
inits_ID <- function(){
	lambda = runif(1, 0.1, 1)
	sigma = runif(1, 1, 5)
	g0 = runif(1,1,10)

	pkj <- (1-exp(-g0*exp(-d2mask/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + 0.0001)%*%t(capt) + log(1-x + 0.00001)%*%t(1-capt))})
	X <- cbind(runif(K, xlim[1], xlim[2]), 
			  runif(K, ylim[1], ylim[2]))
			  
	for(i in 1:K){
		pID <- panimal[ID == i, ]
		if(sum(ID == i) == 1){ 
			mpt <- sample(ncol(panimal), 1, prob = exp(panimal[ID == i, ]))
		}else{
			mpt <- sample(ncol(panimal), 1, prob = exp(colSums(panimal[ID == i, ])))
		}
		X[i,] <- as.numeric(mask[mpt,])
	}
	sigmatoa = runif(1, 0.8, 1)
	list(
        lambda = lambda,
        sigma = sigma,
		sigmatoa = sigmatoa,
		g0 = g0,
		X=X
    )
}

inits_IDSC <- function(){
	lambda = runif(1, 0.1, 1)
	sigma = runif(1, 1, 5)
	g0 = runif(1,1,10)

	pkj <- (1-exp(-g0*exp(-d2mask/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + 0.0001)%*%t(capt) + log(1-x + 0.00001)%*%t(1-capt))})
	X <- cbind(runif(K, xlim[1], xlim[2]), 
			  runif(K, ylim[1], ylim[2]))
			  
	for(i in 1:K){
		pID <- panimal[ID == i, ]
		if(sum(ID == i) == 1){ 
			mpt <- sample(ncol(panimal), 1, prob = exp(panimal[ID == i, ]))
		}else{
			mpt <- sample(ncol(panimal), 1, prob = exp(colSums(panimal[ID == i, ])))
		}
		X[i,] <- as.numeric(mask[mpt,])
	}
	sigmatoa = runif(1, 0.8, 1)
	list(
        lambda = lambda,
        sigma = sigma,
		sigmatoa = sigmatoa,
		g0 = g0,
		X=X,
		psi = rbeta(1,1,1)
    )
}

inits_IDDA <- function(){
	lambda = runif(1, 0.1, 1)
	sigma = runif(1, 1, 5)
	g0 = runif(1,1,10)
	psi = rbeta(1,1,1)

	pkj <- (1-exp(-g0*exp(-d2mask/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x+0.00001)%*%t(capt) + log(1-x+0.00001)%*%t(1-capt))})
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	z = c(rep(NA, K), rep(0, M-K))
	for(i in 1:K){
		pID <- panimal[ID == i, ]
		if(sum(ID == i) == 1){ 
			mpt <- sample(ncol(panimal), 1, prob = exp(panimal[ID == i, ]))
		}else{
			mpt <- sample(ncol(panimal), 1, prob = exp(colSums(panimal[ID == i, ])))
		}
		X[i,] <- as.numeric(mask[mpt,])
	}
	sigmatoa = runif(1, 0.8, 1)
	list(
        lambda = lambda,
        sigma = sigma,
		sigmatoa = sigmatoa,
		g0 = g0,
		X=X,
		psi = psi,
		z = z
    )
}

#########################
## Running ID-ASCR models
#########################

#########################################
## 1) Conditional Likelihood Approach:
#########################################

Rmodel <- nimbleModel(ascr_id_cond, constants_ID, data_ID, inits = inits_ID())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('lambda', 'sigma','sigmatoa', 'g0', 'panimal', 'ESA'))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Make sure it runs...
# Cmcmc$run(10000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out.d <- mcmc(samples[-(1:1000), c("lambda", "sigma", "sigmatoa", "g0", "D", "panimal")])
# plot(out.d)
# hist(out.d[, "D"])
# hist(K/out.d[,'panimal']/area, add = TRUE, col = 'red')

# hist((K/out.d[,'panimal'])/area, freq= FALSE, xlim = c(150,600))
# hist(do.call('c', lapply(out.d[,'panimal'], FUN = function(x){(K+rnbinom(1, size = K, x))/area})), freq=FALSE, add= TRUE, col = 'red')

# ESA <- out.d[,'panimal']*area
# hist(KSA)

## Will need to add the correct error term to density here. e.g. rpois(K)/ESA_hat
t.d <- system.time(
	 out.d <- runMCMC(Cmcmc, niter=50000, nburnin=10000, thin=1, nchains=3,
		inits = list(inits_ID(), inits_ID(), inits_ID()))
)

# tmp <- out.d
# var1 <- K*(1-out.d[[1]][, 'panimal'])
# var2 <- K*(1-out.d[[2]][, 'panimal'])
# var3 <- K*(1-out.d[[3]][, 'panimal'])

# out.d[[1]] <- as.mcmc(cbind(out.d[[1]], "D" = rnorm(nrow(out.d[[1]]), K, sqrt(var1))/out.d[[1]][,"ESA"]/2))
# out.d[[2]] <- as.mcmc(cbind(out.d[[2]], "D" = rnorm(nrow(out.d[[2]]), K, sqrt(var2))/out.d[[2]][,"ESA"]/2))
# out.d[[3]] <- as.mcmc(cbind(out.d[[3]], "D" = rnorm(nrow(out.d[[3]]), K, sqrt(var3))/out.d[[3]][,"ESA"]/2))
# out.d <- as.mcmc.list(out.d)

## Add Poisson variation to K
out.d[[1]] <- as.mcmc(cbind(out.d[[1]], "D" = rpois(nrow(out.d[[1]]), K)/out.d[[1]][,"ESA"]/2))
out.d[[2]] <- as.mcmc(cbind(out.d[[2]], "D" = rpois(nrow(out.d[[2]]), K)/out.d[[2]][,"ESA"]/2))
out.d[[3]] <- as.mcmc(cbind(out.d[[3]], "D" = rpois(nrow(out.d[[3]]), K)/out.d[[3]][,"ESA"]/2))
out.d <- as.mcmc.list(out.d)

# user  	system  elapsed 
# 3820.28    2.78 3825.41
# 1658.08   45.44  440.79 ## C++

# save(out.d, file = "../output/MCMCConditionalDensityCPP.Rda")
# load("../output/MCMCConditionalDensity.Rda")
effectiveSize(out.d)
efficiency.d <- effectiveSize(out.d)/t.d[1]

summary(out.d)

#########################################
## 2) Semi-Complete Likelihood Approach
#########################################

Rmodel <- nimbleModel(ascr_id_SC, constants_IDSC, data_IDSC, inits = inits_IDSC())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('lambda', 'sigma','sigmatoa', 'g0', 'D', 'psi'))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## Make sure it runs...
# Cmcmc$run(10000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out.dsc <- mcmc(samples[-(1:1000), c("lambda", "sigma", "sigmatoa", "g0", "D")])
# plot(out.dsc)
# hist(out.dsc[, "D"])


t.dsc <- system.time(
	 out.dsc <- runMCMC(Cmcmc, niter=50000, nburnin=10000, thin=1, nchains=3,
		inits = list(inits_IDSC(), inits_IDSC(), inits_IDSC()), samplesAsCodaMCMC = TRUE)
)

# user  	system  elapsed 
# 3889.78   3.30    3895.06
# 1684.59   37.88    440.80  ## C++
# save(out.dsc, file = "../output/MCMCFrogSemiCompleteDensityCPP.Rda")
effectiveSize(out.dsc)
efficiency.d <- effectiveSize(out.dsc)/t.dsc[1]

#########################################
## 3)Complete Likelihood Approach w/ data augmentation
#########################################

Rmodel <- nimbleModel(ascr_id_DA, constants_IDDA, data_IDDA, inits = inits_IDDA())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('lambda', 'sigma','sigmatoa', 'g0', 'N', 'D'))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## Make sure it runs...
# Cmcmc$run(10000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out.dda <- mcmc(samples[-(1:1000), c("lambda", "sigma", "sigmatoa", "g0", "D", 'N')])
# plot(out.dda)
# plot(out.dda[,'N'])
# hist(out.dda[, "D"])


t.dda <- system.time(
	 out.dda <- runMCMC(Cmcmc, niter=50000, nburnin=10000, thin=1, nchains=3,
		inits = list(inits_IDDA(), inits_IDDA(), inits_IDDA()), samplesAsCodaMCMC = TRUE)
)

## user    system  elapsed 
## 539.94  0.42    540.48
# save(out.dda, file = "../output/MCMCFrogDataAugmentedDensity.Rda")

effectiveSize(out.dda)/t.dda[1]
summary(out.dda)


hist(out.dda[, "D"])
hist(out.dsc[, "D"], add = TRUE, col = 'red')
hist(out.d[, "D"], add = TRUE, col = 'blue')

MCMCglmm::posterior.mode(out.d[, "D"], 1)
MCMCglmm::posterior.mode(out.dsc[, "D"], 1)
MCMCglmm::posterior.mode(out.dda[, "D"], 1)

effectiveSize(out.d)
effectiveSize(out.dsc)
effectiveSize(out.dda)

par(mfrow = c(2,1))
traceplot(out.d[, 'D'], main = "Conditional Likelihood", xlab = "iters", ylab = "Density", ylim = c(200,700))
traceplot(out.dsc[, 'D'], main = "Semi-Complete Likelihood", xlab = "iters", ylab = "Density", ylim = c(200,700))
traceplot(out.dda[, 'D'], main = "Complete Likelihood w/ Data Augmentation", xlab = "iters", ylab = "Density", ylim = c(200,700))
dev.off()

