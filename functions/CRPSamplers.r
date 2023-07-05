##--------------------------------------------------------------------------
## Nimble Implementation of the CRP
## This is Neal's Algorithm 8 for SCR.
##--------------------------------------------------------------------------
sampler_myCRP <- nimbleFunction(
    name = 'sampler_myCRP',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		m <- control$m
        probs <- numeric(M)
        logProbs <- numeric(M)		
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		keep <- 0
		
		currentValue <- model[[target]]
		n <- sum(model[['ID']] == currentValue)
		if(n == 1){
			keep <- -1
			logProbs[currentValue] <<- model$getLogProb(calcNodes) + log(model[['alpha']]/m)
		}else{
			logProbs[currentValue] <<- model$getLogProb(calcNodes) + log(n-1)
		}
        for(k in 1:M) {
			if(k != currentValue)
			{
				if(model[['z']][k] != 0)
				{
					n <- sum(model[['ID']] == k)
					model[[target]] <<- k
					logProbs[k] <<- model$calculate(calcNodes) + log(n)
				}else{
					logProbs[k] <<- -Inf
				}
			}
		}
		
		a <- m + keep
		for(k in 1:a)
		{
			newVal <- rcat(1, 1-model[['z']])
			model[[target]] <<- newVal
			logProbs[newVal] <<- model$calculate(calcNodesNoSelf) + log(model[['alpha']]/m)			
		}
		
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
		model[[target]] <<- newValue
        if(newValue != currentValue) {
			model[['z']][newValue] <<- 1
			if(keep == -1) model[['z']][currentValue] <<- 0			
            model$calculate(nodesAll)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = nodesAll, logProb = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

##--------------------------------------------------------------------------
## Nimble Implementation of alpha gamma mixture from West 1992
##--------------------------------------------------------------------------
sampler_myAlpha <- nimbleFunction(
    name = 'sampler_myAlpha',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
		shape <- extractControlElement(control, 'shape', 1)		
		rate <- extractControlElement(control, 'rate', 1)		
		nObs <- length(model[['ID']])
    },
    run = function() {
		# This is described in the dirichletprocess package vignette and from West 1992
		# We assume a Gamma prior with parameters alpha_0 (user input).
		# Draw random gamma eta:
		alphai <- model[['alpha']]
		eta <- rbeta(1, alphai + 1, nObs)
		nGroups <- sum(model[['z']])
		
		# Calculate mixing proportions:
		pi1 <- shape + nGroups - 1
		pi2 <- nObs*(rate - log(eta))
		piSelect <- pi1/(pi1 + pi2)

		# Sample from posterior alpha:
		if(runif(1) < piSelect){
			model[[target]] <<- rgamma(1, shape + nGroups, rate - log(eta))
		}else{
			model[[target]] <<- rgamma(1, shape + nGroups - 1, rate - log(eta))
		}
		model$calculate(calcNodes)
		
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list( reset = function() { } )
)

##--------------------------------------------------------------------------
## Sampler for animal activity centres for when they exist and returns fast if not.
##--------------------------------------------------------------------------
sampler_myzs_CRP <- nimbleFunction(
    name = 'sampler_myzs_CRP',
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
## Sampler for animal activity centres for when they exist and returns fast if not.
## This sampler samples z_k = 0 from the conditional distribution of calls.
##--------------------------------------------------------------------------
sampler_myzs_CRP_cond <- nimbleFunction(
    name = 'sampler_myzs_CRP',
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
				lpcurrent <- model[['p.']][nodeIndex]
				model$calculate(calcNodes)
				lpprop <- model[['p.']][nodeIndex]
				logMHR <- lpprop - lpcurrent
				jump <- decide(logMHR)
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