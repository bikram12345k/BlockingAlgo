## This code recreates the results of Figure 1 and 
##  Table 1 from Section 3.3 of the paper.

## Simulation setting
source('bottlenecknbpv2.R')
library('nbpMatching')
library(mixtools)

set.seed(0841)
n = 100

# maxvalues = matrix(NA, 0, 2)
# meanvalues = matrix(NA, 0, 2)

# maxvaluesratio = matrix(NA, 0, 2)
# meanvaluesratio = matrix(NA, 0, 2)


# max5values = matrix(NA, 0, 2)

# distances = matrix(NA, 0, 2)
ITR = 100

est1 = est2 = est3 = est4 = array(0, c(ITR, 5000, 2)) # 4 methods, 2 treatment effect models
ate = matrix(0, ITR, 2) # ATE for 100 iterations and homogeneous or heterogeneous treatment effect
							# models

for(itr in 1:ITR){

		##scenario 1		
		# w = matrix(c(runif(n^2-1000, 0, 2), runif(1000, 4, 5)), n, n)
		# w = w - diag(diag(w))
		# w = (w + t(w))/2

		# x1 = sum(eigen(w)$vectors[,c(1,2)])^2
		# x2 = sum(eigen(w)$vectors[,c(n,n-1,n-2)]/3)
		
		
		##scenario 2
		#library(mixtools)
		#x = rmvnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = rbind(c(0,0), c(2, 2), c(4, 5)), sigma = matrix(rep(1,6),3,2))
		#x1 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
		#x2 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))

		##scenario 3
		d = 6
		Sig = matrix(.3, d, d)
		diag(Sig) = 1
		
		library(mvtnorm)
		dt <- mvtnorm::rmvnorm(n, mean = rep(0,d), sigma = Sig)
		## simulate components
		compidx <- sample(1:3, n, replace=TRUE, prob=c(1/2,3/8,1/8))
		compmeans = matrix(c(0, 2, 5, 0, 2, 5, 0, -2, -6, 0, -2, -6), ncol=3, byrow=TRUE)
		compmeans = rbind(0, 0, compmeans)

		dt = t(compmeans[,compidx])+dt

		dt[,1] = 1*(dt[,1]<.1)
		dt[,2] = 1*(dt[,2]>.1)
		
		
		wrank <- smahal(c(rep(0, n), rep(1, n)), rbind(dt, dt))
		wrank <- sqrt(wrank)
		w <- mahal(dt)
		w <- sqrt(w)
		
		
		w = wrank
		rownames(w) = colnames(w) = 1:dim(w)[1]
		rownames(wrank) = colnames(wrank) = 1:dim(wrank)[1]
		
		
		#w = (outer(x1,x1,'-'))^2 + (outer(x2,x2,'-'))^2
		#w = w#/(2*2)
		#w = sqrt(w)

		rownames(w) = colnames(w) = 1:n

		MAPD = nonbimatch(distancematrix(w))$halves[,c(1,3)]	
		M1 = minmaxnbp(w)

		w_caliper = w
		w_caliper[w_caliper>quantile(as.numeric(w), .90)] = 100*max(w)
		MAPD_caliper = nonbimatch(distancematrix(w_caliper))$halves[,c(1,3)]


		w_caliper = w
		w_caliper[w_caliper>quantile(as.numeric(w), .85)] = 100*max(w)
		MAPD_caliper1 = nonbimatch(distancematrix(w_caliper))$halves[,c(1,3)]

		# outcome model Scenario 2
		# y0 = x1 + x2 + rnorm(n,0,1)
		# y1 = y0 + .25
		# # for scenario 2 y1heter = y0 + abs(x1)
		# y1heter = y0 + sd(y0)*abs(x1)
		
		# outcome model Scenario 3
		y0 = dt[,3]^2*dt[,2] + dt[,4]^2*dt[,1] + dt[,5] + dt[,6] + rnorm(n,0,1)
		y1 = y0 + .5 =
		y1heter = y1 + abs(dt[,4]) 
		 
		 
		ate[itr,] = c(mean(y1 - y0), mean(y1heter - y0))
		
		for(i in 1:5000){
			## Algorithm 1
			trt = runif(n/2, 0, 1)
			trt = 1*(trt<0.5)
			trt = sapply(1:(n/2), function(i) M1[i,1+trt[i]])
			trt1 = rep(0, n)
			names(trt1) = 1:n
			trt1[trt] = 1
			
			bl = rep(NA, n)
			names(bl) = 1:n
			bl[M1[,1]] = 1:(n/2)
			bl[M1[,2]] = 1:(n/2)

			# ATE estimation			
			Y1 = y0*(1-trt1) + y1*trt1
			est1[itr,i,1] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			Y1 = y0*(1-trt1) + y1heter*trt1
			est1[itr,i,2] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			## MAPD
			trt = runif(n/2, 0, 1)
			trt = 1*(trt<0.5)
			trt = sapply(1:(n/2), function(i) MAPD[i,1+trt[i]])
			trt1 = rep(0, n)
			names(trt1) = 1:n
			trt1[trt] = 1
			
			bl = rep(NA, n)
			names(bl) = 1:n
			bl[MAPD[,1]] = 1:(n/2)
			bl[MAPD[,2]] = 1:(n/2)

			# ATE estimation			
			Y1 = y0*(1-trt1) + y1*trt1
			est2[itr,i,1] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			Y1 = y0*(1-trt1) + y1heter*trt1
			est2[itr,i,2] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			## MAPD w/ caliper 10%
			trt = runif(n/2, 0, 1)
			trt = 1*(trt<0.5)
			trt = sapply(1:(n/2), function(i) MAPD_caliper[i,1+trt[i]])
			trt1 = rep(0, n)
			names(trt1) = 1:n
			trt1[trt] = 1
			
			bl = rep(NA, n)
			names(bl) = 1:n
			bl[MAPD_caliper[,1]] = 1:(n/2)
			bl[MAPD_caliper[,2]] = 1:(n/2)

			# ATE estimation			
			Y1 = y0*(1-trt1) + y1*trt1
			est3[itr,i,1] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			Y1 = y0*(1-trt1) + y1heter*trt1
			est3[itr,i,2] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			
			## MAPD w/ caliper 15%
			trt = runif(n/2, 0, 1)
			trt = 1*(trt<0.5)
			trt = sapply(1:(n/2), function(i) MAPD_caliper1[i,1+trt[i]])
			trt1 = rep(0, n)
			names(trt1) = 1:n
			trt1[trt] = 1
			
			bl = rep(NA, n)
			names(bl) = 1:n
			bl[MAPD_caliper1[,1]] = 1:(n/2)
			bl[MAPD_caliper1[,2]] = 1:(n/2)

			# ATE estimation			
			Y1 = y0*(1-trt1) + y1*trt1
			est4[itr,i,1] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
			Y1 = y0*(1-trt1) + y1heter*trt1
			est4[itr,i,2] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
			
		}
	
		
		cat(" - ")
		if(itr %% 10==0)
			cat("\n\n")

	
		## Comparison
		#cat(itr, ': \n')
		# print(c( max(res_sum$halves[,'Distance']), 
				# mean(res_sum$halves[,'Distance']) ) )

		# dist_max = apply(M, 1, function(b) w[b[1], b[2]])

		# print(c(max(dist_max), mean(dist_max)))

		# cat('\n')
		
		# maxvalues = rbind(maxvalues, c(max(res_sum$halves[,'Distance']), 1))
		# maxvalues = rbind(maxvalues, c(max(dist_max), 2))

		# maxvaluesratio = rbind(maxvaluesratio, c(max(res_sum$halves[,'Distance'])/max(dist_max), itr))
		# meanvaluesratio = rbind(meanvaluesratio, c(mean(res_sum$halves[,'Distance'])/mean(dist_max), itr))

		# #print(c(itr, maxvaluesratio[itr,1]*100,  meanvaluesratio[itr,1]*100))
	
		# meanvalues = rbind(meanvalues, c(mean(res_sum$halves[,'Distance']), 1))
		# meanvalues = rbind(meanvalues, c(mean(dist_max), 2))

		# distances = rbind(distances, cbind(res_sum$halves[,'Distance'], 1))
		# distances = rbind(distances, cbind(dist_max, 2))

		# max5values = rbind(max5values, cbind(sort(res_sum$halves[,'Distance'], decreasing=TRUE)[1:2], 1))
		# max5values = rbind(max5values, cbind(sort(dist_max, decreasing=TRUE)[1:2], 2))

}


est1_temp1 = est1_temp2 = est2_temp1 = est2_temp2 =
est3_temp1 = est3_temp2 = est4_temp1 = est4_temp2 = matrix(0, 100, dim(est4)[2])

for(i in 1:dim(est4)[2]){
	est1_temp1[,i] = est1[,i,1] - ate[,1]
	est1_temp2[,i] = est1[,i,2] - ate[,2]
	
	est2_temp1[,i] = est2[,i,1] - ate[,1]
	est2_temp2[,i] = est2[,i,2] - ate[,2]
	
	est3_temp1[,i] = est3[,i,1] - ate[,1]
	est3_temp2[,i] = est3[,i,2] - ate[,2]
	
	est4_temp1[,i] = est4[,i,1] - ate[,1]
	est4_temp2[,i] = est4[,i,2] - ate[,2]
}

ITR1 = itr-1
round(c(mean(apply(est1_temp1[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est2_temp1[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est3_temp1[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est4_temp1[1:ITR1,], 1, function(x) sqrt(mean(x^2))))), 3)


round(c(mean(apply(est1_temp2[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est2_temp2[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est3_temp2[1:ITR1,], 1, function(x) sqrt(mean(x^2)))),
mean(apply(est4_temp2[1:ITR1,], 1, function(x) sqrt(mean(x^2))))), 3)



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


