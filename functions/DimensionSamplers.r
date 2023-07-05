#------------------------------------------------------------------------
## RJMCMC under a fixed parameter space for the purpose of implement in Nimble.
## Using data augmentation indicator variables to do this.
## Simplest verison assumes a discrete uniform prior on N ~ unif(0, M).
#------------------------------------------------------------------------
sampler_myRJMCMC_z <- nimbleFunction(
			name = 'sampler_myRJMCMC_z',
			contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
				zNodes <- target
				zNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		

				loc <- control$ActivityCentre
				delta <- control$delta
				xlim <- control$xlim
				ylim <- control$ylim
				
				unif <- rep(1, delta)
				
				indices <-  as.numeric(gsub(".*?([0-9]+).*", '\\1', zNodesAsScalar))	
				sNodes <- paste0(loc, '[', indices, ', 1:2]')
				calcNodes <- model$getDependencies(c(zNodes, sNodes))				
			},
    run = function() {
		change <- rcat(1,unif)	# Trick where I'm indexing (-delta,..., -1, 1, .. ,delta)
		if(runif(1,0,1) < 0.5) change <- -1*change
		
        currentIndicatorValues <- model[[zNodes]]
		allIndicatorValues <-  model[['z']]	# Ugly hard code quick and dirty.
        currentLogProb <- getLogProb(model, calcNodes)
		NCurrent <- sum(currentIndicatorValues)
		# Can't kill off observed animals or zero animals...
		if(NCurrent + change < 0) return()
		
		## TRUE when kill an observed animal defined by not part of target z.
		killObs <- FALSE
		
		for(k in 1:abs(change))
		{
			if(!killObs){
				## Death
				##------------------------------------------------------------------
				if(change < 0)
				{
					## Careful when you kill you must kill anyone alive which includes 
					## the observed animals! You know that occurs with prob zero so stop 
					## trying if you end up accidentally killing an observed animal for the sake of speed.
					index <- rcat(1, allIndicatorValues)
					if(sum(index == indices) == 0) {
						killObs <- TRUE
					}else{
						# zNodek <- zNodesAsScalar[index]	# Not sure why this notation wasn't working...
						allIndicatorValues[index] <- 0
						# model[[zNodek]] <<- 0
						model[['z']][index] <<- 0
					}
				## Birth
				##------------------------------------------------------------------
				}else{
					index <- rcat(1, 1-allIndicatorValues)
					# zNodek <- zNodesAsScalar[index]
					allIndicatorValues[index] <- 1	
					# model[[zNodek]] <<- 1					## Would be good to get this notation sorted so if someone uses a different letter.
					model[['z']][index] <<- 1

					## Hard coded for uniform priors. Not ideal but a good starting place.
					model[[loc]][index,1] <<- runif(1, xlim[1], xlim[2])
					model[[loc]][index,2] <<- runif(1, ylim[1], ylim[2])
				}
			}
		}

		if(killObs){
			## Use this because our calcNodes don't include observed animals.
			jump <- FALSE
		}else{	
			# Using a uniform proposal distribution 
			proposalLogProb <- calculate(model, calcNodes)
			logMHR <- proposalLogProb - currentLogProb

			jump <- decide(logMHR)
		}
		if(jump)
			copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
		else
			copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
	methods = list(
		reset = function () {}
	)
)

##------------------------------------------------------------------------	
## The point of this sampler is to only sample animals that exist.
## Quick return if the animal doesn't exist for an efficient MCMC algorithm.
##--------------------------------------------------------------------------
sampler_myzs_RJMCMC <- nimbleFunction(
    name = 'sampler_myzs_RJMCMC',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
		scale <- extractControlElement(control, 'scale', 1)		
        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target))		
    },
    run = function() {
		if(model[['z']][nodeIndex] == 0) return()
		newX <- model[[target]] + rnorm(2, 0, scale)
        lpcurrent <- model$getLogProb(calcNodes)
        model[[target]] <<- newX
        lpprop <- model$calculate(calcNodes)
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
			nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
			nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( 
		reset = function() {} 
	)
)

##------------------------------------------------------------------------
## BirthDeathMCMC under a fixed parameter space for the purpose of implement in Nimble.
## Using data augmentation indicator variables to do this.
## Simplest version assumes a discrete uniform prior on N ~ unif(0, M).
## See Stephens (2000) for explanation of method, seems to be a specific form of RJMCMC.
## Birth rate and t0 are both arbitrary. Recommended for t0 = 1 and to choose a reasonable birth rate
## bigger means better mixing and slower. Smaller means faster at the cost of mixing!
##------------------------------------------------------------------------
sampler_myBD_z <- nimbleFunction(
			name = 'sampler_myBD_z',
			contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
				zNodes <- target
				zNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
				calcNodes <- model$getDependencies(zNodes)		

				t0 <- control$t0
				lam <- control$BirthRate
				obsz <- control$ObservedZ
						
				M <- length(obsz)
			},
    run = function() {

		## For internal purposes keep a birth history
		birthHist <- 1-model[['z']]
	
		## Algorithm 3.1 death rate for each component.
		## 2) 
		delta_j <- lam*exp(model[['Hk']])*model[['z']]*(1-obsz)			## Birth Rate * Prob of removing.	

		ti <- 0
		while(ti < t0 & (sum(birthHist) != 0)){
			N <- sum(model[["z"]])
			## 3)
			delta <- sum(delta_j)/N
			## 4)
			trate <- (delta + lam)
			ti <- ti + rexp(1, trate)
			if(ti < t0){
				## 5)
				birth <- rbinom(1, size = 1, prob = lam/trate)
				## 6)
				if(birth == 1){
					index <- rcat(1, birthHist)
					model[['z']][index] <<- 1
					delta_j[index] <- lam*exp(model[['Hk']][index])
				}else{
					index <- rcat(1, delta_j)
					model[['z']][index] <<- 0
					delta_j[index] <- 0
					birthHist[index] <- 0	# Don't birth a dead node, not very independent eh?
				}
			}	
		}
		model$calculate(calcNodes)
		copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
	methods = list(
		reset = function () {}
	)
)
			
##--------------------------------------------------------------------------
## Sampler for activity centres under birth-death sampler.
## The point of this sampler is to only sample animals that exist.
## Quick return if the animal doesn't exist for an efficient MCMC algorithm.
##--------------------------------------------------------------------------
sampler_myzs_BD <- nimbleFunction(
    name = 'sampler_myzs_BD',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
		scale <- extractControlElement(control, 'scale', 1)		
        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target))		
		xlim <- control$xlim
		ylim <- control$ylim
    },
    run = function() {
		if(model[['z']][nodeIndex] == 0) 
		{
				model[[target]][1] <<- runif(1, xlim[1], xlim[2])
				model[[target]][2] <<- runif(1, ylim[1], ylim[2])
				model$calculate(calcNodes)
				jump <- TRUE
		}else{
			newX <- model[[target]] + rnorm(2, 0, scale)
			lpcurrent <- model$getLogProb(calcNodes)
			model[[target]] <<- newX
			lpprop <- model$calculate(calcNodes)
			logMHR <- lpprop - lpcurrent
			jump <- decide(logMHR)
		}
		if(jump) {
			nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
			nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( 
		reset = function() {} 
	)
)

##--------------------------------------------------------------------------
## Stephens 2000 Birth-Death Sampler for SC
##--------------------------------------------------------------------------
sampler_myBDSC_z <- nimbleFunction(
			name = 'sampler_myBDSC_z',
			contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
				zNodes <- target
				zNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
				calcNodes <- model$getDependencies(zNodes)		
				
				t0 <- control$t0
				lam <- control$BirthRate
				M <- control$M
				nobs <- control$ObservedCounts
	
				J <- length(nobs)
			},
    run = function() {

		## For internal purposes keep a birth history
		birthHist <- 1-model[['z']]
		delta_j <- numeric(M)
		Hj <- model[["Hj"]]

		ti <- 0
		while(ti < t0 & (sum(birthHist) != 0)){
			N <- sum(model[["z"]])
			ll.cur <- sum(nobs*log(Hj[1:J])) - sum(Hj[1:J])
			if(N == 0)
			{
				birth <- 1
			}else{
				## Algorithm 3.1 death rate for each component.
				## 2) 
				for(i in 1:M)
				{
					if(model[['z']][i] == 1){
						Hjstar <- Hj[1:J]-model[['Hkj']][i,1:J]
						tmp <- sum(nobs*log(Hjstar)) - sum(Hjstar) - ll.cur
						delta_j[i] <- exp(tmp)*lam/N
					}else{
						delta_j[i] <- 0
					}
				}
				## 3)
				delta <- sum(delta_j)
				## 4)
				trate <- (delta + lam)
				ti <- ti + rexp(1, trate)
				## 5)
				birth <- rbinom(1, size = 1, prob = lam/trate)
			}
			if(ti < t0){
				## 6)
				if(birth == 1){
					index <- rcat(1, birthHist)
					model[['z']][index] <<- 1
					Hj[1:J] <- Hj[1:J] + model[['Hkj']][index,1:J]
				}else{
					index <- rcat(1, delta_j)
					model[['z']][index] <<- 0
					birthHist[index] <- 0	# Don't birth a dead node, not very independent eh?
					Hj[1:J] <- Hj[1:J] - model[['Hkj']][index,1:J]
				}
			}	
		}
		model$calculate(calcNodes)
		copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
	methods = list(
		reset = function () {}
	)
)

##--------------------------------------------------------------------------
## The point of this sampler is to sample both z and s in one step.
## Why do that computation twice? It should make things waaay faster!
##--------------------------------------------------------------------------
sampler_myzs_DA <- nimbleFunction(
    name = 'sampler_myzs_DA',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
		scale <- extractControlElement(control, 'scale', 1)		
        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target[1]))
    },
    run = function() {
        lpcurrent <- model$getLogProb(calcNodes)
		newX <- model[['X']][nodeIndex,] + rnorm(2, 0, scale)
        model[['X']][nodeIndex,] <<- newX
		model[['z']][nodeIndex] <<- 1-model[['z']][nodeIndex]
        lpprop <- model$calculate(calcNodes)
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
			nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
			nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( 
		reset = function() {} 
	)
)
