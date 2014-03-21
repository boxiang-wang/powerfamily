#################################################################################
############ for the use of source ################
#################################################################################
rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")
source("D:/GitHub/powerfamily/U_tool.R")

#################################################################################
############ construct KKT tables ################
#################################################################################

set.seed(1234)
dat = FHTgen(n=100, p=5000, rho=0.8)

system.time(KKTtb(dat, lambda2=0, qv=0.5, nm="n100p5Kr08l20q05"))[3]
system.time(KKTtb(dat, lambda2=1, qv=0.5, nm="n100p5Kr08l21q05"))[3]
system.time(KKTtb(dat, lambda2=0, qv=1, nm="n100p5Kr08l20q1"))[3]
system.time(KKTtb(dat, lambda2=1, qv=1, nm="n100p5Kr08l21q1"))[3]
system.time(KKTtb(dat, lambda2=0, qv=2, nm="n100p5Kr08l20q2"))[3]
system.time(KKTtb(dat, lambda2=1, qv=2, nm="n100p5Kr08l21q2"))[3]


# Summarize the results
 
nmlist = c("n100p5Kr08l20q05", "n100p5Kr08l21q05",
           "n100p5Kr08l20q1", "n100p5Kr08l21q1",
           "n100p5Kr08l20q2", "n100p5Kr08l21q2")

file.loc = "D:\\GitHub\\powerfamily_output\\KKT_check\\"

load(paste(file.loc, nmlist[1], ".rda", sep=""))
m1 = perct.tb
load(paste(file.loc, nmlist[3], ".rda", sep=""))
m3 = perct.tb
load(paste(file.loc, nmlist[5], ".rda", sep=""))
m5 = perct.tb

load(paste(file.loc, nmlist[2], ".rda", sep=""))
m2 = perct.tb
load(paste(file.loc, nmlist[4], ".rda", sep=""))
m4 = perct.tb
load(paste(file.loc, nmlist[6], ".rda", sep=""))
m6 = perct.tb

cbind(rbind(m1, m3, m5), rbind(m2, m4, m6))








#################################################################################
############ others ################
#################################################################################

eps=1e-8
thr=1e-4
lambda2=1
qv=0.5
m.temp = GCDpower(x=dat$x, y=dat$y,
                  lambda2=lambda2, qv=qv, method="power",eps=eps, standardize=F)
KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps, thr=thr)