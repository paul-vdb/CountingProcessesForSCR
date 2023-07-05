simASCR <- function(N = 50, sigma = 0.5, sigma_toa = 0.01, g0 = 1, lambda = 0.5, StudyPeriod = 25, traps, limits, psex = 0.5, nu = 330)
{
    locs <- cbind(x = runif(N, limits[['xlim']][1], limits[['xlim']][2]), 
                  y = runif(N, limits[['ylim']][1], limits[['ylim']][2]))
    J <- nrow(traps)
    obs <- data.frame()
	pop <- data.frame()
	capthist <- toa <- NULL
    ID = 1
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        pkj <- 1-exp(-g0*exp(-d2/(2*sigma^2)))
		p. <- 1-prod(1-pkj)
        ti <- cumsum(rexp(1000, lambda))   #Simulate detection times.
        ti <- ti[ti < StudyPeriod]
        nk <- length(ti)
		keep <- 0
		if(nk != 0) {
			# Now assign those detection times to a trap.
			capt <- do.call("rbind", lapply(1:nk, FUN = function(x){rbinom(J,1,pkj)}))
			keep <- rowSums(capt) != 0
			capt <- capt[keep, ]
			if(sum(keep) > 0) { 
				obs <- rbind(obs, data.frame('t_obs' = ti[keep],'ID' = ID,
														'sex' = rbinom(1, 1, psex),
														'tru_x' = as.numeric(locs[i,1]), 'tru_y' = as.numeric(locs[i,2])))
				capthist <- rbind(capthist, capt)
				
				toa.i <- do.call("rbind", lapply(ti[keep], FUN = function(x){x + sqrt(d2)/nu + rnorm(J,0,sigma_toa)}))
				toa <- rbind(toa, toa.i*capt)
			}
		}
		nk <- sum(keep)
		pop <- rbind(pop, data.frame('n_obs' = nk,'ID' = ID*(nk > 0),
									'tru_x' = locs[i,1], 'tru_y' = locs[i,2]))
		if(nk > 0) 	ID <- ID + 1 							
    }
    list(capt = capthist, toa = toa, obs = obs, pop = pop)
}
