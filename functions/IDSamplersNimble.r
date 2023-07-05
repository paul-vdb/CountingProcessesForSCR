# Small changes here.
# Tried to make it multi-session general but needs some work.
# Seems to be fine for now.
sampler_myBinary <- nimbleFunction(
    name = 'sampler_myBinary',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        nodeIndex <- as.numeric(sub("\\D*(\\d+).*", "\\1", target))
		nv <- extractControlElement(control, 'Noccasion', 1)		
		occ <- extractControlElement(control, 'IDoccasion', 1)
		n_obs <- length(model[['ID']])
		if(nv == 1){ 
			IDMatch <- 1:n_obs
		}else{
			v <- as.numeric(sub(".*,\\D*(\\d+).*", "\\1", target))
			IDMatch <- which(occ == v)
		}
    },
    run = function() {
        if(model[[target]] == 1){
			if(any(model[['ID']][IDMatch] == nodeIndex)) return()
        }
		currentLogProb <- model$getLogProb(calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProb <- model$calculate(calcNodes)
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
        if(jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

######################################################
# Basic sampler for categorical distribution
# I've added the pID here of z as a prior so
# I can do a very quick check on pID[k] = 0
# via calculate(target).
# This generalizes really nicely.
#######################################################
sampler_myCategorical <- nimbleFunction(
    name = 'sampler_myCategorical',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)
    },
    run = function() {
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(k in 1:M) {
            if(k != currentValue) {
                model[[target]] <<- k
                logProbPrior <- model$calculate(target)
                if(logProbPrior == -Inf) {
                    logProbs[k] <<- -Inf
                } else {
                    if(is.nan(logProbPrior)) {
                        logProbs[k] <<- -Inf
                    } else {
                        logProbs[k] <<- logProbPrior + model$calculate(calcNodesNoSelf)
                        if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
                    }
                }
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

######################################################
# Basic sampler for categorical distribution
# sped up for when z is not specifically included.
# Also adding constraints:
# AnimalLink - Observations cannot link to animal index k,
# 	This information may be because we know that animals 1-10 were
#	collared and this observation was not collared.
# cannotlink - Pairwise constraints on the detections themselves that cannot 
# 	be assigned to the same animal.
#######################################################
sampler_mySPIM <- nimbleFunction(
    name = 'sampler_mySPIM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)		
		cannotlink <- extractControlElement(control, 'cannotlink', 'identity')	# n x n matrix where a 1 indicates i cannot link to detection j.
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', targetNodesAsScalar[1]))
		n_grp <- length(targetNodesAsScalar)
    },
    run = function() {
        currentValue <- model[["ID"]][nodeIndex]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(k in 1:M) {
			model[[target]] <<- k
			logProbPrior <- model$calculate(target)
			if(logProbPrior == -Inf) {
				logProbs[k] <<- -Inf
			} else {
				# Add in the pairwise constraint. 1 means cannot link.
				no_link <- sum(cannotlink[model[['ID']] == k, nodeIndex])	# This is the check to see if there are cannot links.
				if(no_link > 0)
				{
					logProbs[k] <<- -Inf
				}else {					
					values(model, targetNodesAsScalar) <<- rep(k, n_grp)
					logProbs[k] <<- model$calculate(calcNodes)
					if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
				}	
			}
		}
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)		
        if(newValue != currentValue) {
			values(model, targetNodesAsScalar) <<- rep(newValue, n_grp)	
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_myCategoricalBernoulli <- nimbleFunction(
    name = 'sampler_myCategoricalBernoulli',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]

		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))

        M <- control$M
		omega <- control$omega
        probs <- numeric(M)
        logProbs <- numeric(M)
    },
    run = function() {
		ID_j <- model[['ID']][omega == omega[nodeIndex]]
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
		
        for(k in 1:M) {
            if(k != currentValue) {
				if(any(ID_j == k)){
					logProbs[k] <<- -Inf
				}else{	
					model[[target]] <<- k
					logProbPrior <- model$calculate(target)
					if(logProbPrior == -Inf) {
						logProbs[k] <<- -Inf
					} else {
						if(is.nan(logProbPrior)) {
							logProbs[k] <<- -Inf
						} else {
							logProbs[k] <<- logProbPrior + model$calculate(calcNodesNoSelf)
							if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
						}
					}
				}
			}	
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

####################################
## Adapted CRP sampler for Fixed dimension
## Detected animals only. New trick for semi-complete sampler.
####################################
sampler_myIDZ <- nimbleFunction(
    name = 'sampler_myIDZ',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		m <- control$m
		
		zNodes <- paste0('z[1:', M, ']')
		Hkv <- 'Hk'
		
        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
		currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodesNoSelf)
		n_currentValue <- sum(model[['ID']] == currentValue)
		keep <- 0
		
		if(n_currentValue == 1)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hkv]][currentValue] 
			keep <- -1
		}else{
			logProbs[currentValue] <<- logProbs[currentValue] + log(1-psi)
		}
				
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[['z']][k] == 0){
					logProbs[k] <<- log(psi) - model[[Hkv]][k]
				}else{
					if(sum(model[['ID']] == k) == 0){
						logProbs[k] <<- -Inf
					}else{
						logProbs[k] <<- log(1-psi)
					}
				}
				model[[target]] <<- k
				logProbs[k] <<- logProbs[k] + model$calculate(calcNodesNoSelf)

				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
					logProbs[k] <<- -Inf
				}
			}
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

####################################
## Adapted CRP sampler for Fixed dimension
## Detected animals only. New trick for semi-complete sampler.
####################################
sampler_myIDZ_cond <- nimbleFunction(
    name = 'sampler_myIDZ_cond',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		m <- control$m
		
		Hk <- 'Hk'
		
        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
		panimal <- model[['panimal']]
		currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodesNoSelf)	## No self as I've put a secret lambda there...
		n_currentValue <- sum(model[['ID']] == currentValue)
		keep <- 0
		
		if(n_currentValue == 1)
		{
			## Change in log prob if we add this one.
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hk]][currentValue] - log(panimal)
			keep <- -1
		}else{
			logProbs[currentValue] <<- logProbs[currentValue] + log(1-psi)
		}
	
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[['z']][k] == 0){
					logProbs[k] <<- log(psi) - model[[Hk]][k] - log(panimal)
				}else{
					logProbs[k] <<- log(1-psi)
				}
				if(logProbs[k] != -Inf){
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodesNoSelf)
				}
				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
					logProbs[k] <<- -Inf
				}
			}
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

####################################
## Adapted CRP sampler for Fixed dimension
## Detected animals only. New trick for semi-complete sampler.
####################################
sampler_myIDZ_SC <- nimbleFunction(
    name = 'sampler_myIDZ_SC',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		m <- control$m
		
		Hk <- 'Hk'
		
        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
		panimal <- model[['panimal']]
		currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodesNoSelf)	## No self as I've put a secret lambda there...
		n_currentValue <- sum(model[['ID']] == currentValue)
		keep <- 0
		
		if(n_currentValue == 1)
		{
			## Change in log prob if we add this one.
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hk]][currentValue]
			keep <- -1
		}else{
			logProbs[currentValue] <<- logProbs[currentValue] + log(1-psi*panimal)
		}
	
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[['z']][k] == 0){
					logProbs[k] <<- log(psi) - model[[Hk]][k]
				}else{
					logProbs[k] <<- log(1-psi*panimal)
				}
				if(logProbs[k] != -Inf){
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodesNoSelf)
				}
				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
					logProbs[k] <<- -Inf
				}
			}
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


####################################
## Adapted CRP sampler for Fixed dimension
## Detected animals only. New trick for semi-complete sampler.
####################################
sampler_myIDZ_SC_8 <- nimbleFunction(
    name = 'sampler_myIDZ_SC_8',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		pchng <- control$pchng
		
		Hk <- 'Hk'
		
        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
		panimal <- model[['panimal']]
		currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodesNoSelf)	## No self as I've put a secret lambda there...
		n_currentValue <- sum(model[['ID']] == currentValue)
		keep <- 0
				
		if(n_currentValue == 1)
		{
			## Change in log prob if we add this one.
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hk]][currentValue]
			keep <- -1
		}else{
			logProbs[currentValue] <<- logProbs[currentValue] + log(1-psi*panimal)
		}
	
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[['z']][k] == 0){
					if(rbinom(1, size = 1, prob = pchng) == 1)
					{
						logProbs[k] <<- log(psi) - model[[Hk]][k]
					}else{
						logProbs[k] <<- -Inf
					}
				}else{
					logProbs[k] <<- log(1-psi*panimal)
				}
				if(logProbs[k] != -Inf){
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodesNoSelf)
				}
				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
					logProbs[k] <<- -Inf
				}
			}
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



####################################
# Generalized for Multi Session:
####################################
sampler_myIDZ_MultiSession <- nimbleFunction(
    name = 'sampler_myIDZ_MultiSession',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		# Find those that are dependent and not dependent for ID
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]

        M <- control$M
		nv <- extractControlElement(control, 'Noccasion', 1)		
		v <- extractControlElement(control, 'occasion', 1)
		if(nv == 1){ 
			zNodes <- paste0('z[1:', M, ']')
			Hkv <- 'Hk'
		}else{
			zNodes <- paste0('z[1:', M,',', v,']')
			Hkv <- paste0('Hk[1:', M, ',', v,']')
		} 
		calcNodesZ <- model$getDependencies(zNodes)
		# Find those that are dependent and not dependent for z
        calcNodesZNoSelf <- model$getDependencies(zNodes, self = FALSE)		
        isStochCalcNodesZNoSelf <- model$isStoch(calcNodesZNoSelf)
        calcNodesZNoSelfDeterm <- calcNodesZNoSelf[!isStochCalcNodesZNoSelf]
        calcNodesZNoSelfStoch <- calcNodesZNoSelf[isStochCalcNodesZNoSelf]

        nodes1 = c(target, zNodes)
        nodes2 = c(calcNodesNoSelfDeterm, calcNodesZNoSelfDeterm)
        nodes3 = c(calcNodesNoSelfStoch, calcNodesZNoSelfStoch)
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		if(n_currentValue == 1)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hkv]][currentValue]
			model[[zNodes]][currentValue] <<- 0
		}
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[[zNodes]][k] == 0){
					logProbs[k] <<- log(psi)-model[[Hkv]][k]				
				}else{
					if(sum(model[['ID']] == k) == 0){
						logProbs[k] <<- -Inf
					}else{
						logProbs[k] <<- log(1-psi)
					}
				}
				# If it is not a zero-inflated value then find the full conditional.
				if(logProbs[k] != -Inf)
				{
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodes)
				}
				
				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
						logProbs[k] <<- -Inf
				}
			}
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
			model[[target]] <<- newValue
            model$calculate(calcNodes)						
			# If I added a new zero animal, or killed an existing one, let's make sure it's updated for z.
			if(model[[zNodes]][newValue] == 0 | model[[zNodes]][currentValue] == 0){
				model[[zNodes]][newValue] <<- 1
				model$calculate(calcNodesZ)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = nodesAll, logProb = TRUE)
			}else{
				model[[zNodes]][currentValue] <<- 1
				nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)				
			}
        } else {
				# Make sure we have the current value switched on!
				model[[zNodes]][currentValue] <<- 1
				model[[target]] <<- currentValue				
				nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
		}
	},
    methods = list( reset = function() { } )
)

# Joint ID Z sampler with constraints.
sampler_myIDZSPIM <- nimbleFunction(
    name = 'sampler_myIDZSPIM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)		
		cannotlink <- extractControlElement(control, 'cannotlink', 'identity')	# n x n matrix where a 1 indicates i cannot link to detection j.
		AnimalLink <- extractControlElement(control, 'AnimalLink', -1)	# A zero represents an impossible animal match.
        if(AnimalLink == -1) AnimalLink <- numeric(length = M, value = 1)
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', targetNodesAsScalar[1]))
		n_grp <- length(targetNodesAsScalar)
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[["ID"]][nodeIndex]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		no_link <- 0
		if(n_currentValue == n_grp)
		{
			model[['z']][currentValue] <<- 0
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi)-model[['Hk']][currentValue]
		}
        for(k in 1:M) {
		if(AnimalLink[k] == 0) {
			logProbs[k] <<- -Inf
		}else {	
				if(k != currentValue){
					if(model[['z']][k] == 1) {
						nk <- sum(model[['ID']] == k)
						if(nk == 0)
						{
							logProbs[k] <<- -Inf
						}else {
							no_link <- sum(cannotlink[model[['ID']] == k, nodeIndex])	# This is the check to see if there are cannot links.
							if(no_link > 0)
							{
								logProbs[k] <<- -Inf
							}else {
								values(model, targetNodesAsScalar) <<- rep(k, n_grp) 					
								logProbs[k] <<- model$calculate(calcNodes) + log(1-psi)
								if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
							}
						}
					}else {
						values(model, targetNodesAsScalar) <<- rep(k, n_grp)
						logProbs[k] <<- model$calculate(calcNodes) + log(psi)-model[['Hk']][k]
						if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
					}
				}
			}
		}
		# Note that logProbs of z=1 and nk=0 is -Inf, or it had better be!y
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
			values(model, targetNodesAsScalar) <<- rep(newValue, n_grp) ##  replace with this			
			if(model[['z']][newValue] == 0){
				model[['z']][newValue] <<- 1
				model$calculate(calcNodesZ)
			}
            model$calculate(calcNodes)	# I've made ID independent of z so this shouldn't double effort.
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
			model[['z']][currentValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)