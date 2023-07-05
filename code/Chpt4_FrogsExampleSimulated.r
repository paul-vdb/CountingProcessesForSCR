########################
# Chapter 4 Examples
# Simulated Frog Acoustic Example
# Latent ID with DPMM
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

library(here)
setwd(here())

# Load Nimble Samplers:
files <- dir("functions")
files <- grep('\\.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0("functions/", x))})
source(paste0("functions", "/", files[6]))

load("data/Frog/stacked-lightfooti.Rdata")
source("SimData.R")

traps.secr <- read.traps(data = data.frame(traps))
mask <- make.mask(traps.secr, buffer = 15, spacing = 0.5)
A <- attr(mask, 'area')
nmask <- nrow(mask)
d2mask <- t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2}))

find_mode_N <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

inits <- function(){
	n <- nrow(capt)
	ID <- numeric(n)
	ID <- 1:n
	z <- numeric(M)
	z[ID] <- 1
	sigma <- runif(1, 3, 5)
	g0 <- runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + .Machine$double.eps)%*%t(capt) + log(1-x + .Machine$double.eps)%*%t(1-capt))})
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))	
	
	for(i in 1:M){
		if(sum(ID == i) == 0) next;
		if(sum(ID == i) == 1) pID <- panimal[which(ID == i), ]
		if(sum(exp(pID), na.rm = TRUE) == 0) next;
		mpt <- sample(ncol(panimal), 1, prob = exp(pID))
		X[i,] <- as.numeric(mask[mpt,])
	}
		
	list(
        sigma = sigma,
		sigmatoa = runif(1, 0.01, 1),
		g0 = g0,
		X = X,
        ID = ID,
		z = z,
		alpha = rgamma(1, shape = 1, rate = 1),
		ESA = 0.02
    )
}

code <- nimbleCode({
	sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
	tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0, 1)
	g0 ~ dunif(0, 20)
	alpha ~ dgamma(shape = 0.1, rate = 0.1)
	psi ~ dbeta(1,1)	
	
	for(k in 1:M) {
		z[k] ~ dbern(0.5)
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		expTime[k, 1:J] <- sqrt(d2[k,1:J])/nu
		pkj[k,1:J] <- (1-exp(-g0*exp(-d2[k,1:J]*tau2)))
		p.[k] <- 1-prod(1-pkj[k,1:J])
	}

	# Trap history model.
	# and unobserved animal ID.
	for(i in 1:n_obs) {
		## Dummy ID distribution to avoid the dCRP...
		ID[i] ~ dID(z[1:M])
		# Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = pkj[ID[i],1:J])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i], 1:J], sd = sigmatoa, y = y[i,1:J])
		# lam[i] <- 1/p.[ID[i]]
	}

	psum <- sum(p.[1:M]*z[1:M])

	ESA <- ESA_C(d2mask[1:nmask, 1:J], sigma, g0, area = A, detfn = 2, nmask)
	pcall <- ESA/area
	## Fast semi-complete addition...
	p0 <- pcall^(-n_obs)
	zero ~ dTrick(p0)
	n ~ dbinom(size = L, prob = pcall*psi)
	
	# Derived Variables.
	K <- sum(z[1:M])
	CD <- psi*L/area
	lam.hat <- n_obs/psum
	D <- CD/lam.hat
})

set.seed(81)
dat <- simASCR(N = 50, sigma = 2.5, sigma_toa = 0.00068, g0 = 5.7, lambda = 0.29, 
	StudyPeriod = 30, traps = traps, 
	limits = list(xlim = range(mask[,1]), ylim = range(mask[,2])))

xlim <- range(mask[,1])
ylim <- range(mask[,2])
area <- nmask*A
toa <- dat$toa
capt <- dat$capt
ID <- dat$obs$ID
IDBen <- ID

# Constants:
M <- 150
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30

constants <- list(
	L = 1000,
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
    n_obs = nrow(capt),
	trials = rep(1, J),
	nu = nu,
	area = area,
	mi = rowSums(capt),
	nmask = nrow(mask),
	d2mask = d2mask,
	A = A
	)

data <- list(
    y = capt,
	toa = toa,
	ID = rep(NA, nrow(capt)),
	z = rep(NA, M),
	zero = 0,
	n=n
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('alpha', 'sigma', 'sigmatoa', 'g0', 'K', 'psi', 'ESA', 'CD', 'D', 'lam.hat'))

conf$removeSamplers('ID')
conf$removeSamplers('z')
# Sampler from Chandler and Royle 2013
conf$addSampler('ID', type = 'myCRP', scalarComponents = TRUE, control = list(M = M, m = 3))

conf$removeSamplers('alpha')
# Sampler from Chandler and Royle 2013
conf$addSampler('alpha', type = 'myAlpha',  control = list(shape = 0.1, rate = 0.1))

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'sampler_myzs_CRP', control = list(scale = 1, xlim=xlim, ylim=ylim))
}

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Make sure it runs...
# Cmcmc$run(10000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)[-(1:5000),]
# out <- as.mcmc(samples)
# plot(out)

out <- runMCMC(Cmcmc, nburnin = 10000, niter = 30000, nchains = 3, 
	inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)
# save(out, file = "output/CRPFrogs_CD_sim.Rda")
# load("../output/CRPFrogs_CD_sim.Rda")

## Semi-Complete ID-ASCR model
ascr_id_SC <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
	g0 ~ dunif(0, 50)
	psi ~ dbeta(1,1)
	
    for(k in 1:K_Obs) {
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

	ESA <- ESA_A(d2mask[1:nmask, 1:J], lambda*Time, sigma, g0, area = A, detfn = 2, nmask)
	panimal <- ESA/area
	p0 <- panimal^(-K)
	zero ~ dTrick(p0)

	## Unobserved Animals:
	# N0 ~ dnegbin(panimal, K)
	K ~ dbinom(size = M, prob = psi*panimal)
	
    # Predicted population size
	D <- (psi*M)/area
})

ID <- IDBen
K <- max(ID)


constants_IDSC <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
	n_obs = nrow(capt),
	trials = rep(1, J),
	nu = nu,
	area = area,
	A = attr(mask, 'area'),
	d2mask = d2mask,
	nmask = nmask,
	ID = IDBen,
	M = M,
	K_Obs = K
	)
##-------------------
data_IDSC <- list(
	zero = 0,
    y = capt,
	toa = toa,
	counts = as.numeric(table(ID)),
	K = K
)
##-------------------

##########################################
## MCMC Initializations
##########################################
inits_IDSC <- function(){
	lambda = runif(1, 0.5, 1)
	sigma = runif(1, 2, 3)
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

Rmodel <- nimbleModel(ascr_id_SC, constants_IDSC, data_IDSC, inits = inits_IDSC())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('lambda', 'sigma','sigmatoa', 'g0', 'D', 'psi'))

conf$removeSamplers('X')
for(i in 1:K){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'),
			type = 'RW_block', silent = TRUE)
}

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out.sc <- runMCMC(Cmcmc, nburnin = 10000, niter = 30000, nchains = 3, 
	inits = list(inits_IDSC(), inits_IDSC(), inits_IDSC()), samplesAsCodaMCMC = TRUE)
# save(out.sc, file = "output/Frogs_ID_sim.Rda")

library(ggplot2)
library(reshape2)

facet_names <- c('D~(Ind~Ha^-1)', 'K', 'sigma~(m)', 'g[0]', 'lambda~(sec^-1)', 'sigma[t]')
names(facet_names) <- c("D", "K", "sigma", "g0", "lambda", "sigmatoa")

out.df <- data.frame(do.call("rbind", out))
out.df$Method = "DPMM"
out.df$lambda <- out.df$lam.hat/30
out.sc.df <- data.frame(do.call("rbind", out.sc))
out.sc.df$Method = "ID-ASCR"
out.sc.df$K = max(IDBen)

tru.vals <- data.frame(variable = c("sigma", "lambda", "g0", "K", "D", "sigmatoa"),
	value = c(2.5, 0.29, 5.7, 14, 50/area, 0.00068))

all.dat <- rbind(out.df[,  c("g0", "sigma", "lambda", "D", "K", "sigmatoa", "Method")], 
	out.sc.df[, c("g0", "sigma", "lambda", "D", "K", "sigmatoa", "Method")])
all.dat.l <- melt(all.dat, id.vars = "Method")
all.dat.l$Method <- factor(all.dat.l$Method)
ggplot(data = all.dat.l, aes(y = value, x = Method)) + 
	facet_wrap(~variable, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) + 
	theme_classic() +
	geom_boxplot(width = 0.075, outlier.alpha = 0, fill = "grey", colour = "black") + 
	geom_violin(alpha = 0.25, adjust = 5, aes(fill = Method),  show.legend = FALSE) + 
	ylab("") + xlab("") + 
	geom_point(data = data.frame(variable = 'K', value = 14, Method = 'ID-ASCR'), shape = 4, size = 3) +
	theme(legend.position = 'bottom') +
	geom_hline(data = tru.vals, aes(yintercept = value), col = 'red') +
	scale_fill_brewer(palette = 1, direction = 1) + 
	theme(axis.text=element_text(size=16),
	axis.title=element_text(size=16),
	legend.title = element_text(size=16), #change legend title font size
	legend.text = element_text(size=16),
	strip.text.x = element_text(size = 16))
ggsave("../output/FrogDPMMResults_SimTOA.png", dpi = 'print', 
	width = 9.75, height = 5.35, units = 'in')

plot(density(do.call('c', out.sc[,'D'])), ylim = c(0,0.01))
lines(density(do.call('c', out[,'D'])), col = 'red')

plot(density(do.call('c', out.sc[,'sigma'])), ylim = c(0,4))
lines(density(do.call('c', out[,'sigma'])), col = 'red')

plot(density(do.call('c', out.sc[,'lambda'])), ylim = c(0,50))
lines(density(do.call('c', out[,'lam.hat'])/30), col = 'red')



ggplot(data = dat$obs, aes(x = tru_x, y = tru_y)) + geom_point() +
	geom_point(data = data.frame(traps), aes(x = x, y = y), col = 'red', shape = 4) + 
	theme_classic() 
