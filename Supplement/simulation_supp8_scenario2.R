rm(list=ls(all=TRUE))

setwd("C:\\Users\\bikramk\\Google Drive\\1- Research and Study\\Research\\Math Stat\\Evidence Factor\\Approx Algo\\code")

library(quickblock)
library(blockTools)
library(optmatch)
library(blockingChallenge)
source('bottlenecknbpv2.R')
library('nbpMatching')
library(mixtools)
dyn.load('derigs_blockingv2.dll')


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



set.seed(1234)

n = 400


ITR = 500-400

fname = "scenario2_k8.txt"
nmethods = 8
#nmethods = 7
k = 8

## The results to be stored.
max_max = matrix(0, nmethods, 0)
max_mean = matrix(0, nmethods, 0)
mean_mean = matrix(0, nmethods, 0)

itr = 1
temp.time = Sys.time()

while(itr<=ITR){
	## Scenario 1
	#w = matrix(c(runif(n^2-1000, 0, 2), runif(1000, 4, 5)), n, n)
	#w = w - diag(diag(w))
	#w = (w + t(w))/2

	## Scenario 2
	x1 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))
	x2 = rnormmix(n, lambda = c(1/2, 3/8, 1/8), mu = c(0, 2, 5))

	w = (outer(x1,x1,'-'))^2 + (outer(x2,x2,'-'))^2
	w = w#/(2*2)
	w = sqrt(w)
			
	rownames(w) = colnames(w) = 1:dim(w)[1]

	## 1
	M1 <- minmaxnbp(w)
	M1vec = apply(M1, 1, FUN=function(m) paste(m, collapse=", "))
	w2 <- outer(M1vec, M1vec, FUN=vneww, w=w)

	M2 <- minmaxnbp(w2)

	rownames(M1) = 1:nrow(M1)
	M2 <- t(apply(M2, 1, function(m) c(M1[m[1],], M1[m[2],])))
	
	M2vec = apply(M2, 1, FUN=function(m) paste(m, collapse=", "))
	w3 <- outer(M2vec, M2vec, FUN=vneww, w=w)

	M3 <- minmaxnbp(w3)
	
	rownames(M2) = 1:nrow(M2)
	design <- t(apply(M3, 1, function(m) c(M2[m[1],], M2[m[2],])))
	
	## 3
	w2 <- outer(M1vec, M1vec, FUN=vneww.1, w=w)
	M2.1 <- minmaxnbp(w2)
	M2.1 <- t(apply(M2.1, 1, function(m) c(M1[m[1],], M1[m[2],])))

	M2vec.1 = apply(M2.1, 1, FUN=function(m) paste(m, collapse=", "))
	w3 <- outer(M2vec.1, M2vec.1, FUN=vneww.1, w=w)

	M3.1 <- minmaxnbp(w3)
	
	rownames(M2.1) = 1:nrow(M2.1)
	design.1 <- t(apply(M3.1, 1, function(m) c(M2.1[m[1],], M2.1[m[2],])))
	
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
	b.k2017 = deprintize(makeblocks)(w, k, maxItr=30)


	# Moore
	## 8
	d = data.frame(id = 1:n, x1 = x1, x2 = x2)
	out <- block(d, n.tr = k, id.vars = c("id"), 
			block.vars = c("x1", "x2"), algorithm="optGreedy", 
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
			if (count == 5) 
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
			if (count == 5) 
				break 
		} else {
			count = 0
			design.1.improve = design.1.improve.temp
		}

	}



	max_max = cbind(max_max, 
					c(max( apply(design, 1, function(m) max(w[m,m])  )),
					max( apply(design.improve, 1, function(m) max(w[m,m])  )),
					max( apply(design.1, 1, function(m) max(w[m,m])  )),
					max( apply(design.1.improve, 1, function(m) max(w[m,m])  )),
					max( apply(b.k2017$strata, 1, function(m) max(w[m,m])  )),
					max(maxd),
					max(maxd.1),
					max( apply(design.qb, 1, function(m) max(w[m,m])  ))
					))


	max_mean = cbind(max_mean, 					
					c(max( apply(design, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					max( apply(design.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					max( apply(design.1, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					max( apply(design.1.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					max( apply(b.k2017$strata, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					max(meand),
					max(meand.1),
					max( apply(design.qb, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  ))
					))

	mean_mean = cbind(mean_mean, 
					c(mean( apply(design, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					mean( apply(design.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					mean( apply(design.1, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					mean( apply(design.1.improve, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					mean( apply(b.k2017$strata, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  )),
					mean(meand),
					mean(meand.1),
					mean( apply(design.qb, 1, function(m) sum(w[m,m])/( length(m)*(length(m)-1) )  ))
					))

					
	writedata = c(itr, max_max[,itr], max_mean[,itr], mean_mean[,itr], 
						mean(table(assg1)), mean(table(assg1.1))) 
	writedata = matrix(writedata, 1)
	
	write.table(writedata, file=fname, col.names = FALSE, row.names = FALSE, append=TRUE)
	
	
	cat(" - ")
	if(itr %% 10==0){
		timetoend = Sys.time() - temp.time
		timetoend = timetoend/10*(ITR-itr)
		temp.time = Sys.time()
		
		cat("\nTime remaining: ", timetoend/(60*60), "hrs. \n\n")
		
	}

	itr = itr + 1
}


