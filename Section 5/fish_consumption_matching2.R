
## Load matching algorithm and code
library(optmatch)
library(approxmatch)
source("nrbalancematch2.R")
assignInNamespace("nrbalancematch", nrbalancematch, "approxmatch")
source("kwaymatching.R")
assignInNamespace("kwaymatching", kwaymatching, "approxmatch")





newvar <- c("age",  "sex", "hispanic", "white", "black", 
			"edu_1", "edu_2", "edu_3", "edu_4", "ipr", "ipr_missing",
			"num_smoke_last30", "smoke_now_everyday", "smoke_now_someday",
			"smoke_100")

rownames(d5) = d5$SEQN
distmat <- multigrp_dist_struc(d5, 'treatment', 
					list(prop=newvar, mahal = newvar), 
					wgts=c(2, 2))

match <- kwaymatching(distmat, 'treatment', .data=d5)

covbalance(d5, 'treatment', match$matches, newvar, details='mean')


d5_match <- d5[rownames(d5) %in% as.vector(match$matches),]

matches <- match$matches
matches <- as.vector(matches)
matches <- data.frame(SEQN = as.numeric(matches), gp = rep(1:nrow(match$matches), 3))

d5_match <- join(d5_match, matches)
