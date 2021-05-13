ITR1=itr-1
tab1 <- array(0, c(ITR1, 3, 6))

for(itr in 1:ITR1){

	tab1[itr,,] = rbind(
colMeans((est1[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est2[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est3[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est4[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est5[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est6[itr,,] - matrix(rep(ate[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2))

}

apply(sqrt(tab1), c(2,3), mean)

tab10 <- array(0, c(ITR, 3, 6))

for(itr in 1:ITR1){

	tab10[itr,,] = rbind(
colMeans((est10[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est20[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est30[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est40[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est50[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2),
colMeans((est60[itr,,] - matrix(rep(ate0[itr,],dim(est1)[2]), dim(est1)[2],byrow=TRUE))^2))

}

apply(sqrt(tab10), c(2,3), mean)


colMeans(ate)

colMeans(est1[1,1:4000,])
ate[1,]







