
d <- read.csv("matched_design2.csv")

names(d)

library(foreign)

d6 <- read.xport('IHGEM_I.xpt')

names(d6)

install.packages('plyr')
library(plyr)

d1<- join(d, d6)

names(d6)

sum(is.na(d1$LBXBGM))

d1$High = d1$fish>=15
d1$Moderate = (d1$fish>=1) & (d1$fish<15)
d1$Low = (d1$fish==0) 

## High consumption to no consumption
senm2(-d1$LBXBGM[d1$Moderate==0],
	d1$High[d1$Moderate==0],
	d1$gp[d1$Moderate==0],
	alternative="less",gamma=4.4)$result$pval

## High consumption to moderate consumption
senm2(-d1$LBXBGM[d1$Low==0],
	d1$High[d1$Low==0],
	d1$gp[d1$Low==0],
	alternative="less",gamma=2.1)$result$pval


senm2(-d1$LBXBGM[d1$High==0],
	d1$Moderate[d1$High==0],
	d1$gp[d1$High==0],
	alternative="greater",gamma=2)$result$pval

d1$treatment.factor = factor(d1$treatment, levels=c(2,3,1),
			labels=c('no', 'between 1 and 5', 'more than 15')) 

par(las=1)
boxplot(log(d1$LBXBGM)~d1$treatment.factor, xlab='Number of Fish Servings in Last 30 Days',
		ylab='Blood methyl-mercury Level, log scale')