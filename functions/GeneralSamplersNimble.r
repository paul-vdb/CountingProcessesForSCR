##--------------------------------------------------------------------------
## Sampler for uniform prior on activity centres. Not in use anymore but
## but an initial attempt to speed up latent ID SCR.
##--------------------------------------------------------------------------
sampler_myX <- nimbleFunction(
    name = 'sampler_myX',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- control$xlim
        ylim <- control$ylim
        calcNodesAll <- model$getDependencies(target)
        calcNodesNoIDs <- model$getDependencies(target)
        # if(model$getDistribution(target) != 'dX') stop('myX sampler error')
        nodeIndex <- as.numeric(gsub('^X\\[([[:digit:]]+), .*$', '\\1', target))
        zNode <- paste0('z[', nodeIndex, ']')
        dNode <- paste0('d2[', nodeIndex, ', 1:', control$J, ']')
        copyNodes <- c(target, dNode)
    },
    run = function() {
        newx <- runif(1, xlim[1], xlim[2])
        newy <- runif(1, ylim[1], ylim[2])
        if(model[[zNode]] == 0) {
            logMHR <- 0
            jump <- decide(logMHR)
            if(jump) {
                model[[target]] <<- c(newx, newy)
                model$calculate(copyNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodes, logProb = FALSE)
            }
            return()
        }
        anyID <- any(model[['ID']] == nodeIndex)
        if(anyID) { lpcurrent <- model$getLogProb(calcNodesAll)
        } else { lpcurrent <- model$getLogProb(calcNodesNoIDs) }
        model[[target]] <<- c(newx, newy)
        if(anyID) { lpprop <- model$calculate(calcNodesAll)
        } else { lpprop <- model$calculate(calcNodesNoIDs) }
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
            if(anyID) {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        } else {
            if(anyID) {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        }
    },
    methods = list( reset = function() {} )
)

##--------------------------------------------------------------------------
# This is my new sampler for multimodal data such as this "label-switching" problem.
# The proposal distribution is a mixture distribution of a uniform and normal.
# The scale is for the normal standard deviation and the
# 'temp' is the temperature, which really is just the mixing proportion of how often to jump.
# If we are running hot, that means mostly uniform, cold is mostly RW.
# Proposal is symmetric so cancels out in MH update.
##--------------------------------------------------------------------------
sampler_myJAM <- nimbleFunction(
    name = 'sampler_myJAM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- extractControlElement(control, 'xlim', c(0,1))		
        ylim <- extractControlElement(control, 'ylim', c(0,1))		
		scale <- extractControlElement(control, 'scale', 1)		
		temp <- extractControlElement(control, 'temp', 0.5)		
		nv <- extractControlElement(control, 'occ', 1)

        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target))
		if(nv == 1){ 
			zNode <- paste0('z[', nodeIndex, ']')
		}else{
			v <- as.numeric(sub(".*,\\D*(\\d+).*", "\\1", target))
			zNode <- paste0('z[', nodeIndex,',', v, ']')
		}		
    },
    run = function() {
		p <- runif(1,0,1)
		if(p < temp)
		{
			# Turn up the heat and try a big step
			newX <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
		}else{
			# Do a local move.
			newX <- model[[target]] + rnorm(2, 0, scale)
		}
		if(model[[zNode]] == 0) {
            logMHR <- 0
            jump <- decide(logMHR)
            if(jump) {
                model[[target]] <<- newX
                model$calculate(calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = FALSE)
            }
            return()
        }
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
    methods = list( reset = function() {} )
)
