####################################
# Basic Example of CRP:
####################################

library(coda)
library(nimble)
library(gtools)
library(here)

n.k.gam <- n.k.const <- NULL
for(i in 1:10000)
{
	lam.gam <- rgamma(1, 5, 1)
	n.k.gam <- c(n.k.gam, rpois(1, lam.gam))
	n.k.const <- c(n.k.const, rpois(1, 5))
}
png("PoisGamExample.png")
hist(n.k.gam, main = "", xlab = expression(n[k]))
hist(n.k.const, add = TRUE, col = 'red')
dev.off()

## CRP
alpha = c(0.5,1,5)
grp <- data.frame(alpha1 = 1, alpha2 = 1, alpha3 = 1)
ngrps <- mean.size <- data.frame()
xlim <- c(-5,5)
ylim <- c(-5,5)
locs <- list(loc1 = data.frame(x = runif(1,xlim), y = runif(1,ylim)),
	loc2 = data.frame(x = runif(1,xlim), y = runif(1,ylim)),
	loc3 = data.frame(x = runif(1,xlim), y = runif(1,ylim)))
for(i in 2:1000)
{
	ngrps <- rbind(ngrps, apply(grp, 2, FUN = function(x){length(unique(x))}))
	mean.size <- rbind(mean.size, apply(grp, 2, FUN = function(x){mean(table(x))}))

	grp.new <- NULL
	for(j in 1:3)
	{
		grp.p <- c(table(grp[,j]), alpha[j])
		new <- sample(length(grp.p), 1, prob=grp.p)
		if(new == length(grp.p)) locs[[j]] <- rbind(locs[[j]],  
					data.frame(x = runif(1,xlim[1], xlim[2]), 
						y = runif(1,ylim[1], ylim[2])))
		grp.new <- c(grp.new, new)
	}
	grp <- rbind(grp, grp.new)
}

png("CRPExample.png")
plot(ngrps[,1], ylim = c(0,40), type = 'l', xlab = "Iteration", ylab = "Number of Groups")
lines(ngrps[,2], col = 'red')
lines(ngrps[,3], col = 'blue')
legend("topright", 
	legend = c(expression(alpha~'='~0.5), expression(alpha~'='~1), expression(alpha~'='~5 )), 
	col = c("black", "red", "blue"),
	lty= 1)
dev.off()

library(ggplot2)
size <- table(grp[,3])
ggplot(data = locs[[3]], aes(x=x, y=y, size=as.numeric(size))) + geom_point() + 
	theme_classic() + labs(size='Number of\nObservations') +
	xlim(xlim) + ylim(ylim)
ggsave("CRPSpatialExample.png", width = 7, height = 6)

hist(size)
	