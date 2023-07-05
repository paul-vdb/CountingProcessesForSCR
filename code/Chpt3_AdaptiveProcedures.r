####################################
# Basic Example of adaptive Procedures:
####################################
library(coda)
library(nimble)
library(here)
setwd(here())

dat <- rpois(100, lambda = 8.5)

interval <- 200
e <- 0.8
AROpt <- 0.44
scale <- 1
lambda <- runif(1,0,10)
ll.cur <- sum(dpois(dat, lambda, log = TRUE))
lambda.res <- NULL
timesAccepted <- 0
timesRan <- 0
scaleHist <- 1
NumberAdapts <- 0
arHist <- 0

for(i in 1:100000)
{
	if((i %% interval) == 0)
	{
		NumberAdapts <- NumberAdapts + 1
		AR <- timesAccepted/timesRan
		arHist <- c(arHist, AR)
		g <- 10/((NumberAdapts + 3)^e)
		adapt <- exp(g*(AR - AROpt))
		scale <- adapt*scale
		scaleHist <- c(scaleHist, scale)
		timesAccepted <- 0
		timesRan <- 0
	}
	lambda.star <- lambda + rnorm(1,0,scale)			# RW new value.
	if(lambda.star < 0){
		ll.new <- -Inf
	}else{
		ll.new <- sum(dpois(dat, lambda.star, log = TRUE))	# Calc new likelihood
	}
	if(exp(ll.new - ll.cur) > runif(1,0,1))
	{
		lambda <- lambda.star	# Define new value.
		ll.cur <- ll.new		# Swap likelihoods
		timesAccepted <- timesAccepted + 1
	}
	lambda.res <- c(lambda.res, lambda)	# Save result.
	timesRan <- timesRan + 1	# Add another run.
}

traceplot(as.mcmc(lambda.res), main = expression(lambda))
png("output/scalehist.png")
	plot(scaleHist, type = 'l', xlab = "Adaptive Step (200 iterations)", ylab = "Scale History")
dev.off()

plot(arHist, type = 'l', xlab = "Adaptive Step (200 iterations)", ylab = "Scale History")


model <- nimbleCode( {
	lambda ~ dunif(0, 20)
	for(k in 1:N){
		n[k] ~ dpois(lambda)
	}
} )


# Run the same model from Burgar et al. Spatial Count on fisher.
#----------------------------------------------------------------------
constants <- list(
    N = length(dat)
	)

data <- list(
		n = dat
	)

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)
Rmodel <- nimbleModel(model, constants, data)
conf <- configureMCMC(Rmodel)
conf$setMonitors('lambda')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(100000)
scaleHist2 <- Cmcmc$samplerFunctions[[1]]$getScaleHistory() 
plot(scaleHist2, type = 'l')

mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out)

#################
# Mixture Model Problems
#################

model <- nimbleCode( {
	lambda1 ~ dunif(0, 20)
	lambda2 ~ dunif(0, 20)
	lambda <- lambda1*0.7 + lambda2*0.3
	
	for(k in 1:N){
		n[k] ~ dpois(lambda)
	}
} )

dat <- rbinom(100, 1, 0.3)
dat <- dat*rpois(100, 15) + rpois(100, 3)*(1-dat)

constants <- list(
    N = length(dat)
	)

data <- list(
		n = dat
	)

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)
Rmodel <- nimbleModel(model, constants, data)
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('lambda1', 'lambda2', 'lambda'))
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(100000)
scaleHist2 <- Cmcmc$samplerFunctions[[2]]$getScaleHistory() 
plot(scaleHist2, type = 'l')

mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out)