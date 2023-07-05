#################################
## ASCR with DPMM
#################################
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)
library(ggplot2)
library(secr)
library(ascr)
library(tidyverse)

# Load all Nimble functions from my source file.
path <- "C:/Users/Paul/Documents/GitHub/PhDNimbleFunctions"
files <- dir(path)
files <- grep('.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0(path, "/", x))})
source(paste0(path, "/", files[1]))

setwd("C:/Users/Paul/Documents/PhD Research/Thesis/Chapter3/Code/DPMM")

## Simulated:
K <- 14
alpha <- 4.5
lams <- rgamma(K, alpha, 0.5)
n <- rpois(K, lams)
n.. <- sum(n)
sd.loc <- 0.6
xlim <- c(-10,10)
ylim <- c(-10,10)

dat <- data.frame()
for(i in 1:K)
{
	if(n[i] > 0){
		mu.tmp <- runif(2, xlim[1], xlim[2])
		dat <- rbind(dat, data.frame(mu1 = mu.tmp[1], mu2 = mu.tmp[2], 
			x1 = rnorm(n[i], mu.tmp[1], sd = sd.loc),  x2 = rnorm(n[i], mu.tmp[2], 
			sd = sd.loc), ID = i))
	}
}

# alpha = 8
# dat <- data.frame()
# for(i in 1:200)
# {
	# if(i > 1){
		# grp.p <- c(table(dat$ID), alpha)
		# new <- sample(length(grp.p), 1, prob=grp.p)
	# }else{
		# new <- 1
		# grp.p <- 1
	# }
	# if(new == length(grp.p)) {
		# dat <- rbind(dat, data.frame(ID = new, mu1 = runif(1, xlim[1], xlim[2]), mu2 = runif(1,ylim[1], ylim[2])))
	# }else{
		# mat <- dat[which(dat$ID == new)[1], ]
		# dat <- rbind(dat, data.frame(ID = new, mu1 = mat[, "mu1"], mu2 = mat[, "mu2"]))
	# }
# }

# save(dat, file = "dataExample.Rda")
load("dataExample.Rda")

ggplot(data = dat, aes(x = x1, y = x2)) + geom_point() +
	geom_point(aes(x = mu1, y = mu2), colour = 'red') + 
	theme_classic() +
	xlab(expression(s[1])) + ylab(expression(s[2])) +
	xlim(xlim + c(-5,5)) + ylim(ylim + c(-5,5)) +
	theme(axis.text=element_text(size=16),
	axis.title=element_text(size=16),
	legend.title = element_text(size=16), #change legend title font size
	legend.text = element_text(size=16))
ggsave("../output/NSExample.png", dpi = 'print', 
	width = 5.75, height = 5.35, units = 'in')

code <- nimbleCode({
	tau ~ dgamma(0.1,0.1)
	alpha ~ dgamma(shape = 0.1, rate = 0.1)

	for(k in 1:M) {
		z[k] ~ dbern(0.5)
		mu[k, 1] ~ dunif(xlim[1], xlim[2])
		mu[k, 2] ~ dunif(ylim[1], ylim[2])
	}

	for(i in 1:n) {
		## Dummy ID distribution to avoid the dCRP...
		ID[i] ~ dID(z[1:M])
		xy[i,1] ~ dnorm(mu[ID[i], 1], tau = tau)
		xy[i,2] ~ dnorm(mu[ID[i], 2], tau = tau)
	}
	
	# Derived Variables.
	K <- sum(z[1:M])
	sigma <- sqrt(1/tau)
})

inits <- function()
{
	sigma <- runif(1,1,3)
	alpha <- runif(1,1,3)
	mu <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
	ID <- do.call('c', lapply(1:n, function(x){sample(M, 1, prob = dnorm(xy[x,1], mu[,1], sigma) * dnorm(xy[x,2], mu[,2], sigma))}))
	z <- rep(0, M)
	z[ID] <- 1
	return(list(tau = 1/sigma^2, sigma=sigma, alpha=alpha, mu=mu, z=z, ID=ID))
}

out.sum <- list()
for( M in c(50, 100, 250, 500, 1000))
{
	# Constants:
	# M <- 100
	n <- nrow(dat)
	xy = cbind(dat$x1, dat$x2)

	constants <- list(
		M = M,
		xlim = xlim,
		ylim = ylim,
		n = n
		)

	data <- list(
		ID = rep(NA, n),
		# ID = dat$ID,
		z = rep(NA, M),
		# z = c(rep(1,20), rep(0, M-20)),
		xy = xy
	)

	Rmodel <- nimbleModel(code, constants, data, inits = inits())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('alpha', 'sigma', 'z', 'K', 'mu', 'ID'))

	conf$removeSamplers('ID')
	conf$removeSamplers('z')
	conf$addSampler('ID', type = 'myCRP', scalarComponents = TRUE, control = list(M = M, m = 3))

	conf$removeSamplers('alpha')
	conf$addSampler('alpha', type = 'myAlpha',  control = list(shape = 0.1, rate = 0.1))

	conf$removeSamplers('mu')
	for(i in 1:M){
		conf$addSampler(target = paste0('mu[', i, ', 1:2]'), 
			type = 'sampler_myzs_CRP', control = list(scale = 1, xlim=xlim, ylim=ylim))
	}

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	## Make sure it runs...
	Cmcmc$run(20000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)[-(1:5000),]
	out <- as.mcmc(samples)
	out.sum[[i]] <- out
}

# save(out.sum, file = "../output/NSMCMCRes.Rda")
load("../output/NSMCMCRes.Rda")

out.id <- out[,grepl('ID', colnames(out))]
avg.nk <- apply(out.id, 1, FUN = function(x){mean(table(x))})
hist(avg.nk,breaks = 50)
abline(v=mean(table(dat$ID)), col='red')

plot(avg.nk,type = 'l')
abline(h=mean(table(dat$ID)), col='red')

summary(out[,c("sigma", "K", "alpha")])
plot(out[,c("sigma", "K", "alpha")])

results <- data.frame()
post.results <- data.frame()
for( M in c(50, 100, 250, 500, 1000))
{
	res <- data.frame(do.call('cbind', summary(out.sum[[M]][,c('K', 'sigma', 'alpha')])))
	res$L <- M
	res$parameter <- rownames(res)
	results <- rbind(results, res)
	out.df <- data.frame(out.sum[[M]][,c('K', 'sigma', 'alpha')])
	out.df$L <- M
	post.results <- rbind(post.results, out.df)
}
results$L <- factor(results$L)
tru.vals <- data.frame(parameter = c('sigma', 'alpha', 'K'), Mean = c(sd.loc, alpha, K))

plot(density(out[,'alpha']))
abline(v = 5, col = 'red')

plot(density(out[,'sigma']))
abline(v = sd.loc, col = 'red')

hist(out[,'K'])
abline(v = K, col = 'red')



plot(dgamma(seq(0.1, 10, by = 0.1), 1, 1))

library(glmmTMB)
glmmTMB(data = dat, x2~1 + (1|ID))


## Test Alpha:
alphai <- 1
shape <- 0.0001
rate <- 0.0001
n <- 125
nGroups <- 20
for(i in 2:1000){

	alpha <- alphai[i-1]
	eta <- rbeta(1, alpha + 1, n)

	# Calculate mixing proportions:
	pi1 <- shape + nGroups - 1
	pi2 <- n*(rate - log(eta))
	piSelect <- pi1/(pi1 + pi2)
	
	if(runif(1) < piSelect){
		alphai <- c(alphai, rgamma(1, shape + nGroups, rate - log(eta)))
	}else{
		alphai <- c(alphai, rgamma(1, shape + nGroups - 1, rate - log(eta)))
	}
}


####################################
## Known ID Version
####################################

## Constants:
M <- max(dat$ID)
n <- nrow(dat)
xy = cbind(dat$x1, dat$x2)

constants <- list(
	M = M,
	xlim = xlim,
	ylim = ylim,
	n = n
	)

data <- list(
	ID = dat$ID,
	z = rep(1, M),
	xy = xy
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('alpha', 'sigma', 'K', 'mu'))

conf$removeSamplers('alpha')
conf$addSampler('alpha', type = 'myAlpha',  control = list(shape = 0.1, rate = 0.1))

conf$removeSamplers('mu')
for(i in 1:M){
	conf$addSampler(target = paste0('mu[', i, ', 1:2]'), 
		type = 'RW_block')
}

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## Make sure it runs...
Cmcmc$run(20000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)[-(1:5000),]
out.id <- as.mcmc(samples)
# save(out.id, file = "../output/NSMCMCResID.Rda")
load("../output/NSMCMCResID.Rda")

out.df <- data.frame(out.id[,c('K', 'sigma', 'alpha')])
out.df$L <- "Known ID"
post.results <- rbind(post.results, out.df)


# res <- data.frame(do.call('cbind', summary(out.id[, c("K", "sigma", "alpha")])))
# res$L <- "Known ID"
# res$parameter <- rownames(res)
# results <- rbind(results, res)

facet_names <- c('K', 'sigma', 'alpha')
names(facet_names) <- c("K", "sigma",  "alpha")

library(reshape2)

tru.vals <- data.frame(variable= c('sigma', 'alpha', 'K'), 
	value = c(sd.loc, alpha, K))

all.dat.l <- melt(post.results, id.vars = "L")
all.dat.l$Scenario <- factor(all.dat.l$L, levels = c("Known ID", "50", "100", "250", "500", "1000"))
ggplot(data = all.dat.l, aes(y = value, x = Scenario)) + 
	facet_wrap(~variable, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) + 
	theme_classic() +
	geom_boxplot(width = 0.075, outlier.alpha = 0, fill = "grey", colour = "black") + 
	geom_violin(alpha = 0.25, adjust = 5, aes(fill = Scenario),  show.legend = FALSE) + 
	ylab("") + 
	geom_hline(data = tru.vals, aes(yintercept = value), col = 'red') + 
	scale_fill_brewer(palette = 1, direction = 1) + 
	theme(legend.position = 'bottom')  +
	theme(axis.text=element_text(size=16),
	axis.title=element_text(size=16),
	legend.title = element_text(size=16), #change legend title font size
	legend.text = element_text(size=16),
	strip.text.x = element_text(size = 16)) +
	theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
	xlab("Scenario")	
ggsave("../output/NSMCMCResults.png", dpi = 'print', width = 9.75, height = 5.35, units = 'in')


ggplot(data = results, aes(y = Mean, x = L)) + 
	facet_wrap(~parameter, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) +
	theme_classic() +
	geom_point() + 
	geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.2) + 
	ylab("") +
	geom_hline(data = tru.vals, aes(yintercept = Mean), col = 'red') +
	theme(legend.position = 'bottom')  +
	theme(axis.text=element_text(size=16),
	axis.title=element_text(size=16),
	legend.title = element_text(size=16), #change legend title font size
	legend.text = element_text(size=16),
	strip.text.x = element_text(size = 16)) +
	theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
	xlab("Scenario")
ggsave("../output/NSMCMCResults.png", dpi = 'print', width = 9.75, height = 5.35, units = 'in')


