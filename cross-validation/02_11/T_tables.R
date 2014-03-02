w = c(1,2,4,17,19,5,7,14,16,11,13)
qvseq = c(0.5, 1, 2, 5, 100)

dat = "arcene"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
load(paste("train_", dat, ".rda", sep=""))
aa2
aaa1 = aa2
load(paste("res_", dat, ".rda", sep=""))
aa
bb1 = aa[w]


dat = "breast"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
load(paste("train_", dat, ".rda", sep=""))
aa2
aaa2 = aa2
load(paste("res_", dat, ".rda", sep=""))
aa
bb2 = aa[w]


dat = "colon"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
load(paste("train_", dat, ".rda", sep=""))
aa2
aaa3 = aa2
load(paste("res_", dat, ".rda", sep=""))
aa
bb3 = aa[w]

dat = "leuk"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
load(paste("train_", dat, ".rda", sep=""))
aa2
aaa4 = aa2
load(paste("res_", dat, ".rda", sep=""))
aa
bb4 = aa[w]

dat = "prostate"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))

load(paste("train_", dat, ".rda", sep=""))
aa2
aaa5 = aa2

load(paste("res_", dat, ".rda", sep=""))
aa
bb5 = aa[w]

rrr = rbind(aaa1, aaa2, aaa3, aaa4, aaa5)
rrr = cbind(c(100,42,62,72,102), c(10000,22283,2000,7128,6033), rrr)
rownames(rrr) = c("Arcene","Breast","Colon","Leukemia","Prostate")
require(xtable)
xtable(rrr, digits=c(0,0,0,2,2,3,3,3))

rr = rbind(bb1, bb2, bb3, bb4, bb5)
rownames(rr) = c("Arcene","Breast","Colon","Leukemia","Prostate")
require(xtable)
xtable(rr, digits=3)

rrr1 = t(rr[,c(2,4,6,8,10)])
matplot(rrr1, type="l")