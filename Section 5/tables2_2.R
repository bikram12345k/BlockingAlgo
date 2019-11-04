library(foreign)
smoke = read.xport('SMQ_I.xpt')

d = read.csv('matched_design2.csv')

matches <- NULL
for(g in 1:max(d$gp)){
	
	temp = d[d$gp==g, c('X','treatment')]
	matches <- rbind(matches, temp$X[order(temp$treatment)])
}

d$edu_5 = 1*(d$edu==5)

library(plyr)
d1 <- join(d, smoke[,c('SEQN', 'SMQ900')])



bage = mean(nulldist_age >= t_age)*100
bsex = mean(nulldist_sex >= t_sex)*100
bedu = mean(nulldist_edu >= t_edu)*100
bipr = mean(nulldist_ipr >= t_ipr)*100
brace = mean(nulldist_race >= t_race)*100
bsmoking = mean(nulldist_smoking_4 >= t_smoking)*100



## balance values
library(approxmatch)

#d$edu_5 = 1*(d$edu==5)

## since d is the matched data set, mean_before and mean_after are the same
tab = covbalance(d1, 'treatment', matches, c('age', 'sex', 'ipr',
							"edu_1", "edu_2", "edu_3", "edu_4", "edu_5",
							"hispanic", "white", "black", "otherrace",
							"smoke_100", 'smoke_now_everyday', 'smoke_now_someday', 
							'num_smoke_last30', 'SMQ900'
							), 'mean')$details$mean_after

tab[2,] = tab[2,]-1
tab[13,] = 1-(tab[13,]-1)
tab[-c(1,3,16),] = tab[-c(1,3,16),]*100
#tab


tab = rbind(tab[1:3,], NA, tab[4:8,], NA, tab[9:12,], NA, tab[13:16,])

rownames(tab) = c('Age', 'Female', 'Income to poverty ratio', 
			'Education',
			'Less than 9th', '9th-11th grade', 'High school graduate', 'Some college', 'College graduate',
			'Race',
			'Hispanic', 'White', 'Black', 'Other race',
			'Smoking status',
			'Smoked at least 100 cigs.', 
			'Smoke now everyday',
			'Smoke now someday',
			'Avg number of cigs. days smoked in past 30 days')
			

tab = cbind(tab, c(bage, bsex, bipr, bedu, rep(NA, 5), brace, rep(NA, 4), bsmoking, rep(NA, 4)))
colnames(tab) = c('More than 15 servings', 'No servings', '1 to 5 servings', 'Balance')
library(xtable)
xtable(tab[,c(1,3,2,4)], digits = 0)



