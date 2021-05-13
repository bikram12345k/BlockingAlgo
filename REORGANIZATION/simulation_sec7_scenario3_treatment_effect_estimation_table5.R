## Simulations of Section 7: Scenario 3 (Table 5)

rm(list=ls(all=TRUE))

#setwd("C:\\Users\\bikramk\\Google Drive\\1- Research and Study\\Research\\Math Stat\\Evidence Factor\\Approx Algo\\code")

library(quickblock)
library(blockTools)
library(optmatch)
library(blockingChallenge)
source('bottlenecknbpv2.R')
library('nbpMatching')
#library(mixtools)
dyn.load('derigs_blockingv2.dll')
source('fixblockTools.R')

# Function for computing 
# rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties 
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance

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


mahal <- function(dt){
	icov = MASS::ginv(cov(dt))
	m = dim(dt)[1]
	out <- matrix(0,m,m)
	for (i in 1:m) 
		out[i,]<-mahalanobis(dt,dt[i,],icov,inverted=T)
	out
}

deprintize<-function(f){
 return(function(...) {capture.output(w<-f(...));return(w);});
}	

neww <- function(x1, x2, w){
	x = strsplit(x1, ", ")[[1]]
	x = c(x, strsplit(x2, ", ")[[1]])
	wx = w[x,x]
	max(wx)
}

vneww <- Vectorize(neww, vectorize.args = c('x1', 'x2'))

neww.1 <- function(x1, x2, w){
	x = strsplit(x1, ", ")[[1]]
	x = c(x, strsplit(x2, ", ")[[1]])
	wx = w[x,x]
	sum(wx)
}

vneww.1 <- Vectorize(neww.1, vectorize.args = c('x1', 'x2'))



set.seed(0841)

n = 200


ITR = 10

fname = "scenario3_k4.txt"
nmethods = 8
#nmethods = 7
k = 4

## The results to be stored.
max_max = matrix(0, nmethods, 0)
max_mean = matrix(0, nmethods, 0)
mean_mean = matrix(0, nmethods, 0)


ate = ate0 = matrix(0, ITR, 3)
est1 = est2 = est3 = est4 = est5 = est6 = array(0, c(ITR, 5000, 3))
est10 = est20 = est30 = est40 = est50 = est60 = array(0, c(ITR, 5000, 3))

itr = 1
while(itr<=ITR){
	## Scenario 1
	#w = matrix(c(runif(n^2-1000, 0, 2), runif(1000, 4, 5)), n, n)
	#w = w - diag(diag(w))
	#w = (w + t(w))/2

	## Scenario 2
	#x1 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
	#x2 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))

	#w = (outer(x1,x1,'-'))^2 + (outer(x2,x2,'-'))^2
	#w = w#/(2*2)
	#w = sqrt(w)
	
	## Scenario 3
	d = 6
	Sig = matrix(.3, d, d)
	diag(Sig) = 1

	#eigen(Sig)$values

	library(mvtnorm)
	dt <- rmvnorm(n, mean = rep(0,d), sigma = Sig)
	## simulate components
	compidx <- sample(1:3, n, replace=TRUE, prob=c(1/2,3/8,1/8))
	compmeans = matrix(c(0, 2, 5, 0, 2, 5, 0, -2, -6, 0, -2, -6), ncol=3, byrow=TRUE)
	compmeans = rbind(0, 0, compmeans)

	dt = t(compmeans[,compidx])+dt

	dt[,1] = 1*(dt[,1]<.1)
	dt[,2] = 1*(dt[,2]>.1)

	## Outcome model
	## 1	
	y1 = dt[,3]^2*dt[,2] + dt[,4]^2*dt[,1] + dt[,5] + dt[,6] + rnorm(n,0,1)
	y2 = y1 + dt[,6] 
	y3 = y1 + dt[,4] 
	y4 = y1 + dt[,3] 
	## 2
	y20 = y1 + dt[,1]*dt[,2]*3
	y30 = y1 + dt[,1]*3
	y40 = y1 + dt[,2]*3
	names(y1) = names(y2) = names(y3) = names(y4) = 1:n
	names(y20) = names(y30) = names(y40) = 1:n

	ate[itr,] = c(mean(y2-y1), mean(y3-y1), mean(y4-y1))
	ate0[itr,] = c(mean(y20-y1), mean(y30-y1), mean(y40-y1))
	
	########################################
	
	wrank <- smahal(c(rep(0, n), rep(1, n)), rbind(dt, dt))
	wrank <- sqrt(wrank)
	w <- mahal(dt)
	w <- sqrt(w)
	
	
	w = wrank
	rownames(w) = colnames(w) = 1:dim(w)[1]
	rownames(wrank) = colnames(wrank) = 1:dim(wrank)[1]
	
	dtrank = apply(dt, 2, rank)
	#View(dtrank)
	cv<-cov(dtrank)
	vuntied<-var(1:n)
	rat<-sqrt(vuntied/diag(cv))
	cv<-diag(rat)%*%cv%*%diag(rat)

	#cv
	dtrank <- data.frame(id = 1:n, dtrank)

	## 1
	M1 <- minmaxnbp(w)
	M1vec = apply(M1, 1, FUN=function(m) paste(m, collapse=", "))
	w2 <- outer(M1vec, M1vec, FUN=vneww, w=w)

	M2 <- minmaxnbp(w2)

	rownames(M1) = 1:nrow(M1)
	design <- t(apply(M2, 1, function(m) c(M1[m[1],], M1[m[2],])))
	
	## 3
	w2 <- outer(M1vec, M1vec, FUN=vneww.1, w=w)
	M2.1 <- minmaxnbp(w2)
	design.1 <- t(apply(M2.1, 1, function(m) c(M1[m[1],], M1[m[2],])))

	w.d = w
	class(w.d) = "distances"

	## 6
	assg1 <- quickblock(w.d, size_constraint=k)
	
	## 7
	assg1.1 <- quickblock(w.d, size_constraint=(k-1))#, break_large_blocks=TRUE)
	#assg1 <- quickblock(w.d, size_constraint=4, break_large_blocks=TRUE)

	maxd = c()
	meand = c()

	st = 0
	while(st<=max(assg1)){
		bl = which(assg1==st)
		nbl = length(bl)
		maxd = c(maxd, max(w[bl,bl]))
		meand = c(meand, sum(w[bl,bl])/( nbl*(nbl-1) ))
		st=st+1
	}


	maxd.1 = c()
	meand.1 = c()

	st = 0
	while(st<=max(assg1.1)){
		bl = which(assg1.1==st)
		nbl = length(bl)
		maxd.1 = c(maxd.1, max(w[bl,bl]))
		meand.1 = c(meand.1, sum(w[bl,bl])/( nbl*(nbl-1) ))
		st=st+1
	}


	#library(optmatch)
	b.k2017 = deprintize(makeblocks)(wrank, k, maxItr=30)


	# Moore
	## 8
	#d = data.frame(id = 1:n, x1 = x1, x2 = x2)
	colnames(cv) = rownames(cv) = names(dtrank)[-1]
	
	out <- block1(dtrank, n.tr = k, id.vars = c("id"), 
			block.vars = names(dtrank)[-1], algorithm="optGreedy", 
			vcov.data = cv,
			distance = "mahalanobis", level.two = FALSE, verbose = FALSE)
	 
	design.qb = out$blocks[[1]][,-(k+1)]
	design.qb = t(apply(design.qb, 1, as.character))

	## Improvements



	## Improve design
	m = nrow(design)
	d = rep(0, m)
	s = 0

	## 2
	count = 0
	design.improve = design
	while(1){
		c0 = max( apply(design.improve, 1, function(m) max(w[m,m])  ))
		## Find the furthest unit in each block
		badx <- c()
		for(i in 1:m){
			b = design.improve[i,]
			b.temp = sapply(b, function(x) {
									b_x = setdiff(b, x)
									max(w[b_x, b_x])})
			#temp = apply(w[b,b], 1, max)
			temp = b[b.temp == min(b.temp)]
			badx <- c(badx, sample(temp)[1])
		}

		design.reduced = t(apply(design.improve, 1, function(x) setdiff(x, badx)))
		## Create new distance matarix
		w.improve <- matrix(NA, m, m)
		for(i in 1:m){
			x = badx[i]
			w.improve[i,] = sapply(1:m, function(j) {
								tempb = c(x, design.reduced[j,])
								max(w[tempb, tempb])
							})
		}

		w.improve = round(w.improve*(10^(6)))

		res = .Fortran('calllbap', n = as.integer(m), c = as.integer(w.improve), spalte = as.integer(d), z = as.integer(s))

		badx = badx[res$spalte]
		design.improve.temp  = t(sapply(1:m, function(i) c(badx[i], design.reduced[i,])))
		
		c1 = max( apply(design.improve.temp, 1, function(m) max(w[m,m])  ))
		#print(res$spalte)
		if(c0==c1){
			count = count + 1
			if (count == 10) 
				break 
		} else {
			count = 0
			design.improve = design.improve.temp
		}

	}

	## 4
	count = 0
	design.1.improve = design.1
	while(1){
		c0 = max( apply(design.1.improve, 1, function(m) mean(w[m,m])  ))	
		## Find the furthest unit in each block
		badx <- c()
		for(i in 1:m){
			b = design.1.improve[i,]
			b.temp = sapply(b, function(x) {
									b_x = setdiff(b, x)
									mean(w[b_x, b_x])})
			#temp = apply(w[b,b], 1, max)
			temp = b[b.temp == min(b.temp)]
			badx <- c(badx, sample(temp)[1])
		}

		design.1.reduced = t(apply(design.1.improve, 1, function(x) setdiff(x, badx)))
		## Create new distance matarix
		w.improve <- matrix(NA, m, m)
		for(i in 1:m){
			x = badx[i]
			w.improve[i,] = sapply(1:m, function(j) {
								tempb = c(x, design.1.reduced[j,])
								sum(w[tempb, tempb])
							})
		}

		w.improve = round(w.improve*(10^(6)))

		res = .Fortran('calllbap', n = as.integer(m), c = as.integer(w.improve), spalte = as.integer(d), z = as.integer(s))

		badx = badx[res$spalte]
		design.1.improve.temp = t(sapply(1:m, function(i) c(badx[i], design.1.reduced[i,])))
		
		c1 = max( apply(design.1.improve.temp, 1, function(m) mean(w[m,m])  ))
		#print(res$spalte)
		if(c0==c1){
			count = count + 1
			if (count == 10) 
				break 
		} else {
			count = 0
			design.1.improve = design.1.improve.temp
		}

	}



	estfn <- function(trt, designt, y2, y3, y4){
			# 1
		trt1 = apply(trt, 1, function(x) which(x==1))
		Yo1 = y1[sapply(1:length(trt1), function(k) designt[k,trt1[k]])]
		# 2
		trt2 = apply(trt, 1, function(x) which(x==2))
		Yo2 = y2[sapply(1:length(trt2), function(k) designt[k,trt2[k]])]
		# 3
		trt3 = apply(trt, 1, function(x) which(x==3))
		Yo3 = y3[sapply(1:length(trt3), function(k) designt[k,trt3[k]])]
		# 4
		trt4 = apply(trt, 1, function(x) which(x==4))
		Yo4 = y4[sapply(1:length(trt4), function(k) designt[k,trt4[k]])]
		
		c(mean(Yo2), mean(Yo3), mean(Yo4))-mean(Yo1)
	}

	## The treatment effect estimation loop
	for(i in 1:5000){
		trt = t(sapply(1:(n/4), function(x) sample(1:4)))
		# given design
		est1[itr,i,] = estfn(trt, design, y2, y3, y4)
		est2[itr,i,] = estfn(trt, design.improve, y2, y3, y4)
		est3[itr,i,] = estfn(trt, design.1, y2, y3, y4)
		est4[itr,i,] = estfn(trt, design.1.improve, y2, y3, y4)
		est5[itr,i,] = estfn(trt, b.k2017$strata, y2, y3, y4)
		est6[itr,i,] = estfn(trt, design.qb, y2, y3, y4)
		
		est10[itr,i,] = estfn(trt, design, y20, y30, y40)
		est20[itr,i,] = estfn(trt, design.improve, y20, y30, y40)
		est30[itr,i,] = estfn(trt, design.1, y20, y30, y40)
		est40[itr,i,] = estfn(trt, design.1.improve, y20, y30, y40)
		est50[itr,i,] = estfn(trt, b.k2017$strata, y20, y30, y40)
		est60[itr,i,] = estfn(trt, design.qb, y20, y30, y40)
		
		#est2[itr,i] = estfn(trt, )
		#est3[itr,i] = estfn(trt, )
		#est4[itr,i] = estfn(trt, )
		
		
	}
	
	
	

	# max_max = cbind(max_max, 
					# c(max( apply(design, 1, function(m) max(w[m,m])  )),
					# max( apply(design.improve, 1, function(m) max(w[m,m])  )),
					# max( apply(design.1, 1, function(m) max(w[m,m])  )),
					# max( apply(design.1.improve, 1, function(m) max(w[m,m])  )),
					# max( apply(b.k2017$strata, 1, function(m) max(w[m,m])  )),
					# max(maxd),
					# max(maxd.1),
					# max( apply(design.qb, 1, function(m) max(w[m,m])  ))
					# ))


	# max_mean = cbind(max_mean, 					
					# c(max( apply(design, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# max( apply(design.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# max( apply(design.1, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# max( apply(design.1.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# max( apply(b.k2017$strata, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# max(meand),
					# max(meand.1),
					# max( apply(design.qb, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  ))
					# ))

	# mean_mean = cbind(mean_mean, 
					# c(mean( apply(design, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# mean( apply(design.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# mean( apply(design.1, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# mean( apply(design.1.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# mean( apply(b.k2017$strata, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					# mean(meand),
					# mean(meand.1),
					# mean( apply(design.qb, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  ))
					# ))

					
	# writedata = c(itr, max_max[,itr], max_mean[,itr], mean_mean[,itr], 
						# mean(table(assg1)), mean(table(assg1.1))) 
	# writedata = matrix(writedata, 1)
	
	# write.table(writedata, file=fname, col.names = FALSE, row.names = FALSE, append=TRUE)
	
	
	cat(" - ")
	if(itr %% 10==0)
		cat("\n\n")

	itr = itr + 1
}


######## Calculation for complete randomization

set.seed(0841)

n = 200


ITR = 10000

fname = "scenario3_k4.txt"
nmethods = 8
#nmethods = 7
k = 4

## The results to be stored.
max_max = matrix(0, nmethods, 0)
max_mean = matrix(0, nmethods, 0)
mean_mean = matrix(0, nmethods, 0)


ate = ate0 = matrix(0, ITR, 3)
est = est0 = matrix(0, ITR, 3)

itr = 1
while(itr<=ITR){
	## Scenario 1
	#w = matrix(c(runif(n^2-1000, 0, 2), runif(1000, 4, 5)), n, n)
	#w = w - diag(diag(w))
	#w = (w + t(w))/2

	## Scenario 2
	#x1 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
	#x2 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))

	#w = (outer(x1,x1,'-'))^2 + (outer(x2,x2,'-'))^2
	#w = w#/(2*2)
	#w = sqrt(w)
	
	## Scenario 3
	d = 6
	Sig = matrix(.3, d, d)
	diag(Sig) = 1

	#eigen(Sig)$values

	library(mvtnorm)
	dt <- rmvnorm(n, mean = rep(0,d), sigma = Sig)
	## simulate components
	compidx <- sample(1:3, n, replace=TRUE, prob=c(1/2,3/8,1/8))
	compmeans = matrix(c(0, 2, 5, 0, 2, 5, 0, -2, -6, 0, -2, -6), ncol=3, byrow=TRUE)
	compmeans = rbind(0, 0, compmeans)

	dt = t(compmeans[,compidx])+dt

	dt[,1] = 1*(dt[,1]<.1)
	dt[,2] = 1*(dt[,2]>.1)

	## Outcome model
	## 1	
	y1 = dt[,3]^2*dt[,2] + dt[,4]^2*dt[,1] + dt[,5] + dt[,6] + rnorm(n,0,1)
	y2 = y1 + dt[,6] 
	y3 = y1 + dt[,4] 
	y4 = y1 + dt[,3] 
	## 2
	y20 = y1 + dt[,1]*dt[,2]*3
	y30 = y1 + dt[,1]*3
	y40 = y1 + dt[,2]*3
	names(y1) = names(y2) = names(y3) = names(y4) = 1:n
	names(y20) = names(y30) = names(y40) = 1:n

	ate[itr,] = c(mean(y2-y1), mean(y3-y1), mean(y4-y1))
	ate0[itr,] = c(mean(y20-y1), mean(y30-y1), mean(y40-y1))
	

	## Randomize
	trt = sample(rep(1:4,each=50))

	Y1 = y1[trt==1]
	Y2 = y2[trt==2]
	Y3 = y3[trt==3]
	Y4 = y4[trt==4]

	Yo1 = y1[trt==1]
	Yo2 = y20[trt==2]
	Yo3 = y30[trt==3]
	Yo4 = y40[trt==4]

	est[itr,] = c(mean(Y2), mean(Y3), mean(Y4))-mean(Y1)
	est0[itr,] = c(mean(Yo2), mean(Yo3), mean(Yo4))-mean(Yo1)
	
	itr = itr + 1
}

round(sapply(1:3, function(i) sqrt( mean( (est[,i] - ate[,i] )^2 ) ) ), 3)
round(sapply(1:3, function(i) sqrt( mean( (est0[,i] - ate0[,i] )^2 ) ) ), 3)

