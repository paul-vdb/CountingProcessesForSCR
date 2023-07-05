##--------------------------------------------------------------------------
## Fast computation of ESA for semi-complete likelihood and conditional!!
## Assuming Bernoulli capture history. In particular this is for acoustics stevenson 2015
##--------------------------------------------------------------------------
RESA_C <- function(d2mask, sigma, g0, area, detfn, nmask) {

		ESA <- 0
		for(i in 1:nmask)
		{
			if(detfn == 1) # Half-normal
			{
				p.mask <- g0*exp(-d2mask[i,]/(2*sigma^2))
			}else{	# Hazard half-normal
				p.mask <- (1-exp(-g0*exp(-d2mask[i,]/(2*sigma^2))))
			}
			p.detect <- 1-prod(1-p.mask)
			ESA <- ESA + p.detect*area	
		}
   		return(ESA)
}

##--------------------------------------------------------------------------
## Fast computation of ESA for semi-complete likelihood and conditional!!
## This one is for a repeated caller at a fixed location from stevenson 2021
##--------------------------------------------------------------------------
RESA_A <- function(d2mask, lambda, sigma, g0, area, detfn, nmask) {

		ESA <- 0
		for(i in 1:nmask)
		{
			if(detfn == 1) # Half-normal
			{
				p.mask <- g0*exp(-d2mask[i,]/(2*sigma^2))
			}else{	# Hazard half-normal
				p.mask <- (1-exp(-g0*exp(-d2mask[i,]/(2*sigma^2))))
			}
			p.detect <- 1-prod(1-p.mask)
			
			ESA <- ESA + (1-exp(-p.detect*lambda))*area	# Assumes lambda is total expected calls.
		}
   		return(ESA)
}

##--------------------------------------------------------------------------
## Fast computation of ESA for semi-complete likelihood and conditional!!
## Assuming Bernoulli capture history. In particular this is for acoustics stevenson 2015
##--------------------------------------------------------------------------
RESA_Cam <- function(d2mask, lambda, sigma, area, nmask) {

		ESA <- 0
		for(i in 1:nmask)
		{
			Hk <- sum(lambda*exp(-d2mask[i,]/(2*sigma^2)))
			ESA <- ESA + (1-exp(-Hk))*area	
		}
   		return(ESA)
}


##--------------------------------------------------------------------------
## Effective sample area for bernoulli detections
##--------------------------------------------------------------------------
ESA_C <- nimbleRcall(function(d2mask = double(2), sigma=double(0), 
								g0=double(0), area = double(0), detfn = double(0), 
								nmask = double(0)){}, 
		Rfun = 'RESA_C',
		returnType = double(0)
		)
##--------------------------------------------------------------------------
## Effective Sample area for Poisson model
##--------------------------------------------------------------------------	
ESA_A <- nimbleRcall(function(d2mask = double(2), lambda = double(0), sigma=double(0), 
								g0=double(0), area = double(0), detfn = double(0), 
								nmask = double(0)){}, 
		Rfun = 'RESA_A',
		returnType = double(0)
		)
		
##--------------------------------------------------------------------------
## Effective Sample area for Poisson model using camera traps.
##--------------------------------------------------------------------------			
ESA_AC <- nimbleRcall(function(d2mask = double(2), lambda = double(0), sigma=double(0), 
								area = double(0), nmask = double(0)){}, 
		Rfun = 'RESA_Cam',
		returnType = double(0)
		)		


library(Rcpp)
library(RcppArmadillo)
sourceCpp("functions/helperCPP.cpp")

ESA_C_CPP <- nimbleRcall(function(d2mask = double(2), sigma=double(0), 
								g0=double(0), area = double(0), detfn = integer(0)){}, 
		Rfun = 'CPP_ESA_C',
		returnType = double(0)
		)
ESA_A_CPP <- nimbleRcall(function(d2mask = double(2), lambda = double(0), sigma=double(0), 
								g0=double(0), area = double(0), detfn = integer(0)){}, 
		Rfun = 'CPP_ESA_A',
		returnType = double(0)
		)
ESA_AC_CPP <- nimbleRcall(function(d2mask = double(2), lambda = double(0), sigma=double(0), 
								area = double(0)){}, 
		Rfun = 'CPP_ESA_Cam',
		returnType = double(0)
		)