
# Effect of sampling error on analysis of
# Within- and across-species correlations
##########################################
# R. Dakin, November 2017

# Suppose we have 3 species,
# each has within-species correlation between two traits
# (e.g., size, speed)
# but species average trait value is the same for all species
# (i.e., no species differences)

# set up the populations:

npop <- 10000
set.seed(101)
sp1 <- data.frame(size=rep(10,npop)+rnorm(npop), species='sp1', stringsAsFactors=F)
sp2 <- data.frame(size=rep(10,npop)+rnorm(npop), species='sp2', stringsAsFactors=F)
sp3 <- data.frame(size=rep(10,npop)+rnorm(npop), species='sp3', stringsAsFactors=F)
sp1$speed <- sp1$size*2 + rnorm(npop,0,1)
sp2$speed <- sp2$size*2 + rnorm(npop,0,1)
sp3$speed <- sp3$size*2 + rnorm(npop,0,1)

head(sp1)

# function to normalize a variable on a scale from 0-1 (e.g., for plotting):
mynorm <- function(x){
	temp <- x - min(x) # set min to 0
	temp <- temp/max(temp) # set max to 1
	return(temp)
}

# plot (some of) the population-level data
dev.new(width=7,height=2.5)
par(mfrow=c(1,3), bty='l', las=1, mar=c(4,4,0.25,0.2))
plot(speed ~ size, sp1[1:100,], pch=16, cex=mynorm(size)+0.25, xlim=c(6,12), ylim=c(12,26))
abline(h=20, v=10, lty=3)
plot(speed ~ size, sp2[1:100,], pch=17, cex=mynorm(size)+0.25, xlim=c(6,12), ylim=c(12,26))
abline(h=20, v=10, lty=3)
plot(speed ~ size, sp3[1:100,], pch=15, cex=mynorm(size)+0.25, xlim=c(6,12), ylim=c(12,26))
abline(h=20, v=10, lty=3)

summary(sp1) # mean size = 10; mean speed = 20 for each species
summary(sp2)
summary(sp3)

# Suppose we draw a small sample from each
# take the species avg. both traits
# and test correlation

nreps <- 5000
ss1 <- 3
results <- data.frame(within.cor=rep(NA, nreps), w.p=NA, between.cor=NA, b.p=NA)
set.seed(202)
for(i in 1:nreps){
	s1 <- sp1[sample(npop, ss1),]
	s2 <- sp2[sample(npop, ss1),]
	s3 <- sp3[sample(npop, ss1),]
	samp <- rbind(s1,s2,s3)
	smeans <- aggregate(samp[,c('size','speed')], by=list(species=samp$species), 'mean')
	names(smeans)[2:3] <- c('sp.size','sp.speed')
	samp <- merge(samp, smeans, by='species')
	samp$rel.size <- samp$size - samp$sp.size
	results$within.cor[i] <- cor.test(samp$speed, samp$rel.size)$estimate
	results$w.p[i] <- cor.test(samp$speed, samp$rel.size)$p.value
	results$between.cor[i] <- cor.test(smeans$sp.speed, smeans$sp.size)$estimate
	results$b.p[i] <- cor.test(smeans$sp.speed, smeans$sp.size)$p.value
	print(i)
}
head(results)

sum(results$b.p<0.05)/length(results$b.p) # around 16%, more than 3x the expected rate of false positives, all because of the within-species correlations (note: the results also includes the analysis of the within-species correlation although we don't look at it here.)

# What if we increase the sample sizes drawn?

nreps <- 5000
ss1 <- 10
set.seed(303)
for(i in 1:nreps){
	s1 <- sp1[sample(npop, ss1),]
	s2 <- sp2[sample(npop, ss1),]
	s3 <- sp3[sample(npop, ss1),]
	samp <- rbind(s1,s2,s3)
	smeans <- aggregate(samp[,c('size','speed')], by=list(species=samp$species), 'mean')
	names(smeans)[2:3] <- c('sp.size','sp.speed')
	samp <- merge(samp, smeans, by='species')
	samp$rel.size <- samp$size - samp$sp.size
	results$within.cor10[i] <- cor.test(samp$speed, samp$rel.size)$estimate
	results$w.p10[i] <- cor.test(samp$speed, samp$rel.size)$p.value
	results$between.cor10[i] <- cor.test(smeans$sp.speed, smeans$sp.size)$estimate
	results$b.p10[i] <- cor.test(smeans$sp.speed, smeans$sp.size)$p.value
	print(i)
}
head(results)

sum(results$b.p10<0.05)/length(results$b.p10) # 16%! this does not solve the problem of inflated false positives. Even increasing it to 100 or even 1000 individuals per species doesn't help.

# Plot the distribution of p-values for the species-level correlation:
dev.new(width=5, height=2.5)
par(mfrow=c(1,2), las=1, mar=c(4,4,2,0.5))
hist(results$b.p, main='Test between sp. corr,\n n = 3 per species', cex.main=0.75, breaks=20, xlim=c(0,1), ylim=c(0,1000), xlab='p-value', ylab='frequency', xaxt='n', yaxt='n')
axis(1, at=c(0, 0.5,1))
axis(2, at=c(0,250,500))
abline(v=0.05, lty=3, col='red')
hist(results$b.p10, main='Test between sp. corr,\n n = 10 per species', cex.main=0.75, breaks=20, xlim=c(0,1), ylim=c(0,1000), xlab='p-value', ylab='frequency', xaxt='n', yaxt='n')
axis(1, at=c(0, 0.5,1))
axis(2, at=c(0,250,500))
abline(v=0.05, lty=3, col='red')

# What if we break the (individual-level) link between speed and size, by drawing one small sample for speed and one small sample for size.

nreps <- 5000
ss1 <- 3
ss2 <- 3
results2 <- data.frame(within.cor=rep(NA, nreps), w.p=NA, between.cor=NA, b.p=NA)
set.seed(404)
for(i in 1:nreps){
	s1 <- sp1[sample(npop, ss1),c('size','speed','species')]
	s2 <- sp2[sample(npop, ss1),c('size','speed','species')]
	s3 <- sp3[sample(npop, ss1),c('size','speed','species')]
	s4 <- sp1[sample(npop, ss2),c('size','species')]
	s5 <- sp2[sample(npop, ss2),c('size','species')]
	s6 <- sp3[sample(npop, ss2),c('size','species')]
	
	samp1 <- rbind(s1,s2,s3)
	smeans <- aggregate(samp1[,c('speed')], by=list(species=samp1$species), 'mean')
	
	samp2 <- rbind(s4,s5,s6)
	smeans2 <- aggregate(samp2[,c('size')], by=list(species=samp2$species), 'mean')
	samp2 <- merge(samp1, smeans2, by='species')
	names(samp2)[4] <- 'sp.size'
	
	samp2$rel.size <- samp2$size - samp2$sp.size
	
	results2$within.cor[i] <- cor.test(samp2$speed, samp2$rel.size)$estimate
	results2$w.p[i] <- cor.test(samp2$speed, samp2$rel.size)$p.value
	results2$between.cor[i] <- cor.test(smeans$x, smeans2$x)$estimate
	results2$b.p[i] <- cor.test(smeans$x, smeans2$x)$p.value
	print(i)
}

sum(results2$b.p<0.05)/length(results2$b.p) # approx. 0.05 as expected

set.seed(505)
ss1 <- 10
ss2 <- 10
for(i in 1:nreps){
	s1 <- sp1[sample(npop, ss1),c('size','speed','species')]
	s2 <- sp2[sample(npop, ss1),c('size','speed','species')]
	s3 <- sp3[sample(npop, ss1),c('size','speed','species')]
	s4 <- sp1[sample(npop, ss2),c('size','species')]
	s5 <- sp2[sample(npop, ss2),c('size','species')]
	s6 <- sp3[sample(npop, ss2),c('size','species')]
	
	samp1 <- rbind(s1,s2,s3)
	smeans <- aggregate(samp1[,c('speed')], by=list(species=samp1$species), 'mean')
	
	samp2 <- rbind(s4,s5,s6)
	smeans2 <- aggregate(samp2[,c('size')], by=list(species=samp2$species), 'mean')
	samp2 <- merge(samp1, smeans2, by='species')
	names(samp2)[4] <- 'sp.size'
	
	samp2$rel.size <- samp2$size - samp2$sp.size
	
	results2$within.cor10[i] <- cor.test(samp2$speed, samp2$rel.size)$estimate
	results2$w.p10[i] <- cor.test(samp2$speed, samp2$rel.size)$p.value
	results2$between.cor10[i] <- cor.test(smeans$x, smeans2$x)$estimate
	results2$b.p10[i] <- cor.test(smeans$x, smeans2$x)$p.value
	print(i)
}

sum(results2$b.p10<0.05)/length(results2$b.p10) # again as expected

dev.new(width=5, height=2.5)
par(mfrow=c(1,2), las=1, mar=c(4,4,2,0.5))
hist(results2$b.p, main='Test between sp. corr,\n n = 3 per species', cex.main=0.75, breaks=20, xlim=c(0,1), ylim=c(0,1000), xlab='p-value', ylab='frequency', xaxt='n', yaxt='n')
axis(1, at=c(0, 0.5,1))
axis(2, at=c(0,250,500))
abline(v=0.05, lty=3, col='red')
hist(results2$b.p10, main='Test between sp. corr,\n n = 10 per species', cex.main=0.75, breaks=20, xlim=c(0,1), ylim=c(0,1000), xlab='p-value', ylab='frequency', xaxt='n', yaxt='n')
axis(1, at=c(0, 0.5,1))
axis(2, at=c(0,250,500))
abline(v=0.05, lty=3, col='red')












