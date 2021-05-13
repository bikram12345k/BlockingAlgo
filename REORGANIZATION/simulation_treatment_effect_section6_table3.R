## Table 3 of Section 6

source('bottlenecknbpv2.R')
#install.packages('nbpMatching')
library('nbpMatching')


smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,,drop=FALSE]
    Xt<-X[z==1,,drop=FALSE]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }


X1 = c(0, 4, 14, 10, 4, 5, 6, 6, 4, 10, 7, 4, 1, 9, 4, 0, 11, 6, 0, 6, 11, 7, 0, 12)
X2 = c(0, 2, 14, 5, 5, 4, 8, 6, 6, 8, 6, 3, 1, 10, 10, 0, 10, 6, 8, 5, 10, 6, 0, 11)

YP0 = c(7, 1, 10, 10, 6, 2, 7, 5, 6, 6, 3, 8, 0, 11, 10, 0, 8, 6, 7, 1, 8, 5, 1, 5)
YM0 = c(2, 0, 13, 0, 0, 0, 1, 2, 0, 0, 2, 1, 0, 1, 13, 0, 0, 4, 0, 1, 5, 7, 0, 0)

#predict(lm(YP~X1))
YP = YP0
YM = YP-1

X = c(X1 + X2)/2
n = length(X)

set.seed(814)
## Simulate many 
idx = which(X<=7)
samp = sample(idx, 60, replace=TRUE)
samp = c(samp, sample(setdiff(1:n, idx), 40, replace=TRUE))

load("C:\\Users\\bikra\\Google Drive\\1- Research and Study\\Research\\Math Stat\\Evidence Factor\\Approx Algo\\code\\Repository\\Supplement\\samp.randomsample")

#X1 = X[samp]

#samp = 1:24

n = length(samp)
#w = outer(X1, X1, FUN = function(x,y) abs(x-y))
X1samp = X1[samp]+rnorm(n, 0, .01)
X2samp = X2[samp]+rnorm(n, 0, .01)

w = smahal(c(rep(0,n), rep(1,n)), 
		rbind(cbind(X1samp, X2samp), cbind(X1samp, X2samp)))

rownames(w) = colnames(w) = 1:n
names(samp) = 1:n
#dim(w)


M1 <- minmaxnbp(w)

#M1
MAPD = nonbimatch(distancematrix(w))$halves[,c(1,3)]
MAPD = as.matrix(MAPD)

#apply(M1, 1, function(b) w[b[1], b[2]])

w_caliper = w
w_caliper[w_caliper>.30] = 100
MAPD_caliper = nonbimatch(distancematrix(w_caliper))$halves[,c(1,3)]


w_caliper = w
w_caliper[w_caliper>.25] = 100
MAPD_caliper1 = nonbimatch(distancematrix(w_caliper))$halves[,c(1,3)]


MAPD = as.matrix(MAPD)
MAPD_caliper = as.matrix(MAPD_caliper)
MAPD_caliper1 = as.matrix(MAPD_caliper1)



#w_caliper = w
#w_caliper[w_caliper>.25] = 100
#MAPD_caliper1 = nonbimatch(distancematrix(w_caliper))$halves[,c(1,3)]

cbind(sort(apply(M1, 1, function(b) w[b[1], b[2]])),
		sort(apply(MAPD, 1, function(b) w[b[1], b[2]])))

boxplot( c(apply(M1, 1, function(b) w[b[1], b[2]]),
		apply(MAPD, 1, function(b) w[b[1], b[2]]),
		apply(MAPD_caliper, 1, function(b) w[b[1], b[2]]),
		apply(MAPD_caliper1, 1, function(b) w[b[1], b[2]])) ~ rep(1:4,each=n/2) )

###############################################
te = c(0, .5, .75)#, 1, 1.25)
rejectionCR = matrix(NA, 0, 2)
estCR = matrix(0, 10000, length(te))
pvalCR.0 = pvalCR.1 = NULL

ate = NULL


## Complete randomization
for(teid in 1:length(te)){
	if(teid==1) YM = YP - .25 #te[teid]
	if(teid==2) YM = pmax(0, YP - X1/15)
	if(teid==3) YM = pmax(0, YP - sqrt(X1)/3)
	#if(teid==4) YM = max(0, YP - pmin(X1, 6))
	
	ate = c(ate, mean(YP-YM))

	pvaCR.0 = pvalCR.1 = c()

	for(i in 1:10000){
		trtunits = sample(1:100)[1:(100/2)]
		trt = rep(0, 100)
		trt[trtunits] = 1
		Y1 = YP[samp]*(1-trt) + YM[samp]*trt

		#install.packages('ri2')
		#library(ri2)

		table_2.2 <- data.frame(d = trt, y = Y1)
		declaration <- declare_ra(N = n, m = n/2)

		#out <- conduct_ri(y ~ d,declaration = declaration,assignment = "d",
		#		sharp_hypothesis = 0,data = table_2.2,p='upper')


		#pvalCR.0 = c(pvalCR.0, summary(out)[1,3])

		table_2.2$R1 = lm(Y1~X1samp+X2samp)$residuals

		#out <- conduct_ri(R1 ~ d,declaration = declaration,assignment = "d",
		#		sharp_hypothesis = 0,data = table_2.2,p='upper')


		#pvalCR.1 = c(pvalCR.1, summary(out)[1,3])

		#pval1 = c(pval1, summary(lm(Y1~trt1+X1))$coefficients[2, 4])
		estCR[i,teid] = mean(Y1[trt==1]) - mean(Y1[trt==0])
	}

	#rejectionCR = rbind(rejectionCR, c(mean(pvalCR.0<.05), mean(pvalCR.1<.05)))
	
}


###############################################
te = c(0, .5, .75)#, 1, 1.25)
rejection = matrix(NA, 0, 4)

est1 = est2 = est3 = est4 = matrix(0, 10000, length(te))

for(teid in 1:length(te)){
	if(teid==1) YM = YP - .25 #te[teid]
	if(teid==2) YM = pmax(0, YP - X1/15)
	if(teid==3) YM = pmax(0, YP - sqrt(X1)/3)
	#if(teid==4) YM = max(0, YP - pmin(X1, 6))
	
	## do many randomization

	pval1.0 = pval1.1 = c()

	for(i in 1:10000){
		trt = runif(n/2, 0, 1)
		trt = 1*(trt<0.5)
		trt = sapply(1:(n/2), function(i) M1[i,1+trt[i]])
		trt1 = rep(0, n)
		names(trt1) = 1:n
		trt1[trt] = 1
		Y1 = YP[samp]*(1-trt1) + YM[samp]*trt1
		bl = rep(NA, n)
		names(bl) = 1:n
		bl[M1[,1]] = 1:(n/2)
		bl[M1[,2]] = 1:(n/2)

		# pval1.0 = c(pval1.0, 
			# as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		# R1 = lm(Y1~X1samp+X2samp)$residuals

		# pval1.1 = c(pval1.1, 
			# as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))

		#pval1 = c(pval1, summary(lm(Y1~trt1+X1))$coefficients[2, 4])
		est1[i,teid] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
	}

	## 
	pval2.0 = pval2.1 = c()

	for(i in 1:10000){

		trt = runif(n/2, 0, 1)
		trt = 1*(trt<0.5)
		trt = sapply(1:(n/2), function(i) MAPD[i,1+trt[i]])
		trt1 = rep(0, n)
		names(trt1) = 1:n
		trt1[trt] = 1
		Y1 = YP[samp]*(1-trt1) + YM[samp]*trt1

		bl = rep(NA, n)
		names(bl) = 1:n
		bl[MAPD[,1]] = 1:(n/2)
		bl[MAPD[,2]] = 1:(n/2)


		# pval2.0 = c(pval2.0, 
			# as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		# R1 = lm(Y1~X1samp+X2samp)$residuals

		# pval2.1 = c(pval2.1, 
			# as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))
		#pval2.1 = c(pval2, summary(lm(Y1~trt1+X1))$coefficients[2,4])
		est2[i,teid] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
	}


	for(i in 1:10000){

		trt = runif(n/2, 0, 1)
		trt = 1*(trt<0.5)
		trt = sapply(1:(n/2), function(i) MAPD_caliper[i,1+trt[i]])
		trt1 = rep(0, n)
		names(trt1) = 1:n
		trt1[trt] = 1
		Y1 = YP[samp]*(1-trt1) + YM[samp]*trt1

		bl = rep(NA, n)
		names(bl) = 1:n
		bl[MAPD_caliper[,1]] = 1:(n/2)
		bl[MAPD_caliper[,2]] = 1:(n/2)


		# pval2.0 = c(pval2.0, 
			# as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		# R1 = lm(Y1~X1samp+X2samp)$residuals

		# pval2.1 = c(pval2.1, 
			# as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))
		#pval2.1 = c(pval2, summary(lm(Y1~trt1+X1))$coefficients[2,4])
		est3[i,teid] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
	}
	
	
	for(i in 1:10000){

		trt = runif(n/2, 0, 1)
		trt = 1*(trt<0.5)
		trt = sapply(1:(n/2), function(i) MAPD_caliper1[i,1+trt[i]])
		trt1 = rep(0, n)
		names(trt1) = 1:n
		trt1[trt] = 1
		Y1 = YP[samp]*(1-trt1) + YM[samp]*trt1

		bl = rep(NA, n)
		names(bl) = 1:n
		bl[MAPD_caliper1[,1]] = 1:(n/2)
		bl[MAPD_caliper1[,2]] = 1:(n/2)


		# pval2.0 = c(pval2.0, 
			# as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		# R1 = lm(Y1~X1samp+X2samp)$residuals

		# pval2.1 = c(pval2.1, 
			# as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))
		#pval2.1 = c(pval2, summary(lm(Y1~trt1+X1))$coefficients[2,4])
		est4[i,teid] = mean(Y1[trt1==1]) - mean(Y1[trt1==0])
	}
	
	rejection = rbind(rejection, c(mean(pval1.0<0.05), mean(pval1.1<0.05), 
					mean(pval2.0<0.05), mean(pval2.1<0.05)))
	#boxplot(c(pval1, pval2)~c(rep(1,1000),rep(2,1000)))

}

# rownames(rejection) = te
# rejection


ate = c()
for(teid in 1:length(te)){
	if(teid==1) YM = YP - .25 #te[teid]
	if(teid==2) YM = pmax(0, YP - X1/15)
	if(teid==3) YM = pmax(0, YP - sqrt(X1)/3)

	ate = c(ate, mean(YP[samp] - YM[samp]))
}

ate

tab <- matrix(0, length(te), 4)
for(teid in 1:length(te)){
tab[teid,] = c(sqrt(mean( (est1[,teid] - ate[teid])^2 )),
	sqrt(mean( (est2[,teid] - ate[teid])^2 )),
	sqrt(mean( (est3[,teid] - ate[teid])^2 )),
	sqrt(mean( (est4[,teid] - ate[teid])^2 )))
}


tab = cbind(tab,  c(sqrt(mean( (estCR[,1] - ate[1])^2 )),
	sqrt(mean( (estCR[,2] - ate[2])^2 )),
	sqrt(mean( (estCR[,3] - ate[3])^2 ))))


colnames(tab) = c('Algo1', 'MAPD', "MAPD_caliper", "MAPD_caliper_repeated", "CR")
rownames(tab) = paste0("Model ", 1:length(te))

tab

library(xtable)
xtable((cbind(ate, tab[,c(1,3,2,5)])), digits=2)







