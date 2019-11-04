## This code recreates the results of Figure 1 and 
##  Table 1 from Section 3.3 of the paper.

## Simulation setting
source('bottlenecknbpv2.R')
library('nbpMatching')

set.seed(0841)
n = 100

maxvalues = matrix(NA, 0, 2)
meanvalues = matrix(NA, 0, 2)

maxvaluesratio = matrix(NA, 0, 2)
meanvaluesratio = matrix(NA, 0, 2)


max5values = matrix(NA, 0, 2)

distances = matrix(NA, 0, 2)

for(itr in 1:100){

		##scenario 1		
		#w = matrix(c(runif(n^2-1000, 0, 2), runif(1000, 4, 5)), n, n)
		#w = w - diag(diag(w))
		#w = (w + t(w))/2


		##scenario 2
		#library(mixtools)
		#x = rmvnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = rbind(c(0,0), c(2, 2), c(4, 5)), sigma = matrix(rep(1,6),3,2))
		x1 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
		x2 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
		
		w = (outer(x1,x1,'-'))^2 + (outer(x2,x2,'-'))^2
		w = w#/(2*2)
		w = sqrt(w)

		rownames(w) = colnames(w) = 1:n

		res_sum = nonbimatch(distancematrix(w))
		M = minmaxnbp(w)

		## Comparison
		cat(itr, ': \n')
		print(c( max(res_sum$halves[,'Distance']), 
				mean(res_sum$halves[,'Distance']) ) )

		dist_max = apply(M, 1, function(b) w[b[1], b[2]])

		print(c(max(dist_max), mean(dist_max)))

		cat('\n')
		
		maxvalues = rbind(maxvalues, c(max(res_sum$halves[,'Distance']), 1))
		maxvalues = rbind(maxvalues, c(max(dist_max), 2))

		maxvaluesratio = rbind(maxvaluesratio, c(max(res_sum$halves[,'Distance'])/max(dist_max), itr))
		meanvaluesratio = rbind(meanvaluesratio, c(mean(res_sum$halves[,'Distance'])/mean(dist_max), itr))

		#print(c(itr, maxvaluesratio[itr,1]*100,  meanvaluesratio[itr,1]*100))
	
		meanvalues = rbind(meanvalues, c(mean(res_sum$halves[,'Distance']), 1))
		meanvalues = rbind(meanvalues, c(mean(dist_max), 2))

		distances = rbind(distances, cbind(res_sum$halves[,'Distance'], 1))
		distances = rbind(distances, cbind(dist_max, 2))

		max5values = rbind(max5values, cbind(sort(res_sum$halves[,'Distance'], decreasing=TRUE)[1:2], 1))
		max5values = rbind(max5values, cbind(sort(dist_max, decreasing=TRUE)[1:2], 2))

}




maxvalues1 = data.frame(maxvalues)
names(maxvalues1) = c('values', 'type')

meanvalues1 = data.frame(meanvalues)
names(meanvalues1) = c('values', 'type')

meanvalues1$type = meanvalues1$type + 2


##scenario 1
#boxplot(values~type, rbind(meanvalues1, maxvalues1))

max5values.sc1 = max5values
aggregate(max5values.sc1[,1], by = list(max5values.sc1[,2]), quantile, c(.5, .75, .95, .99, 1))

maxvaluesratio.sc1 = maxvaluesratio
mean(maxvaluesratio.sc1[,1]*100)

distances.sc1 = distances
boxplot(distances.sc1[,1]~distances.sc1[,2])



distances.sc1.p$Algorithm = factor(distances.sc1.p$Algorithm, levels=c(1,2, 3), labels=c("'optimal' matching", "Algorithm 1", "Unif"))

distances.sc1.1 = distances.sc1
distances.sc1.1 = rbind(distances.sc1.1, cbind(runif(20000, 0, 2), 3))

boxplot(distances.sc1.1[,1]~distances.sc1.1[,2])



##scenario 2
max5values.sc2 = max5values
aggregate(max5values.sc2[,1], by = list(max5values.sc2[,2]), quantile, c(.5, .75, .95, .99, 1))

maxvaluesratio.sc2 = maxvaluesratio
mean(maxvaluesratio.sc2[,1]*100)


#distances.sc2 = distances
distances.sc2.p$Algorithm = factor(distances.sc2.p$Algorithm, levels=c(1,2, 3), labels=c("'optimal' matching", "Algorithm 1", "Chisq"))

distances.sc2.1 = distances.sc2
distances.sc2.1 = rbind(distances.sc2.1, cbind(sqrt(rchisq(20000, 1)), 3))

boxplot(distances.sc2.1[,1]~distances.sc2.1[,2])



qqplot(distances[distances[,2]==2,1], sqrt(rchisq(1000, 2)), xlim=c(0,4), ylim=c(0,4))
#abline(a=0,b=1)
par(new=TRUE)
qqplot(distances[distances[,2]==1,1], sqrt(rchisq(1000, 2)), xlim=c(0,4), ylim=c(0,4), col='grey')
abline(a=0,b=1)




################################# Plots and table

tab = rbind(cbind(Scenario='Scenario 1', aggregate(max5values.sc1[,1], by = list(max5values.sc1[,2]), quantile, c(.5, .75, .95, .99, 1))),
	cbind(Scenario='Scenario 2', aggregate(max5values.sc2[,1], by = list(max5values.sc2[,2]), quantile, c(.5, .75, .95, .99, 1))))

tab = cbind(tab, pct.max=c(mean(1/maxvaluesratio.sc1[,1]*100),100, mean(1/maxvaluesratio.sc2[,1]*100), 100))

tab = data.frame(tab[,-c(2)])
tab

library(xtable)
xtable(tab)



distances.sc1.p = data.frame(distances.sc1)
names(distances.sc1.p) = c("Distances", "Algorithm")
distances.sc1.p$Algorithm = factor(distances.sc1.p$Algorithm, levels=c(1,2), labels=c("MAPD", "Algorithm 1"))

distances.sc1.p.1 = rbind(distances.sc1.p, cbind(Distances=abs(rnorm(20000, 0, 1)), Algorithm=rep(3,20000)))
boxplot(distances.sc1.p.1[,1]~distances.sc1.p.1[,2])


qqplot(distances.sc2.p[distances.sc2.p$Algorithm=="Algorithm 1",1], abs(rnorm(20000, 0, 1)), xlim=c(0,4.5), ylim=c(0,4.5),
		xlab='', ylab='')
par(new=TRUE)
qqplot(distances.sc2.p[distances.sc2.p$Algorithm=="'optimal' matching",1], abs(rnorm(20000, 0, 1)), col='grey', xlim=c(0,4.5), ylim=c(0,4.5),
		xlab='', ylab='')
abline(a=0,b=1, lty=2)


distances.sc1.p = data.frame(distances.sc1)
names(distances.sc1.p) = c("Distances", "Algorithm")
distances.sc1.p$Algorithm = factor(distances.sc1.p$Algorithm, levels=c(1,2), labels=c("MAPD", "Algorithm 1"))

distances.sc2.p = data.frame(distances.sc2)
names(distances.sc2.p) = c("Distances", "Algorithm")
distances.sc2.p$Algorithm = factor(distances.sc2.p$Algorithm, levels=c(1,2), labels=c("MAPD", "Algorithm 1"))

par(mfrow=c(1,2))
par(oma=c(0,1,0,0), las=1, mar=c(5.1,4.1,2.1,1.1))
boxplot(distances.sc1.p[,1]~distances.sc1.p[,2], ylab = 'Difference', ylim=c(0,4.4))
mtext("Scenario 1", side=1, line=2.75)
boxplot(distances.sc2.p[,1]~distances.sc2.p[,2], ylab = 'Difference', ylim=c(0,4.4))
mtext("Scenario 2", side=1, line=2.75)

#mtext("Scenario 1\t\tScenario 2", outer=TRUE, side=1,line=0)


