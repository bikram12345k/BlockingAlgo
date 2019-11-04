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

## Simulate many 
idx = which(X<=7)
samp = sample(idx, 60, replace=TRUE)
samp = c(samp, sample(setdiff(1:n, idx), 40, replace=TRUE))

#X1 = X[samp]

samp = 1:24

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

#apply(M1, 1, function(b) w[b[1], b[2]])

cbind(sort(apply(M1, 1, function(b) w[b[1], b[2]])),
		sort(apply(MAPD, 1, function(b) w[b[1], b[2]])))




###############################################
te = c(0, .5, .75, 1, 1.25)
rejection = matrix(NA, 0, 4)

for(teid in 1:length(te)){
	YM = YP - te[teid]
	## do many randomization

	pval1.0 = pval1.1 = c()

	for(i in 1:4000){
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

		pval1.0 = c(pval1.0, 
			as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		R1 = lm(Y1~X1samp+X2samp)$residuals

		pval1.1 = c(pval1.1, 
			as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))

		#pval1 = c(pval1, summary(lm(Y1~trt1+X1))$coefficients[2, 4])
	}

	## 
	pval2.0 = pval2.1 = c()

	for(i in 1:4000){

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


		pval2.0 = c(pval2.0, 
			as.numeric(pairwisePermutationSymmetry(x=Y1, g=trt1,b=as.factor(bl))[1,3]))

		R1 = lm(Y1~X1samp+X2samp)$residuals

		pval2.1 = c(pval2.1, 
			as.numeric(pairwisePermutationSymmetry(x=R1, g=trt1,b=as.factor(bl))[1,3]))
		#pval2.1 = c(pval2, summary(lm(Y1~trt1+X1))$coefficients[2,4])
	}

	rejection = rbind(rejection, c(mean(pval1.0<0.05), mean(pval1.1<0.05), 
					mean(pval2.0<0.05), mean(pval2.1<0.05)))
	#boxplot(c(pval1, pval2)~c(rep(1,1000),rep(2,1000)))

}

rownames(rejection) = te
rejection




