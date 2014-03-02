# source("D:/GitHub/powerfamily/U_tool.R")
# datapath = "D:/GitHub/powerfamily/data/"
# dyn.load("D:/GitHub/powerfamily/O_hsvmlassoNET.dll")

source("/home/wang3660/Research/PF/U_tool_server.R")
datapath = "/home/wang3660/Research/PF/data/"
dyn.load("/home/wang3660/Research/PF/crossvalidation/O_hsvmlassoNET.so")

dat = "arcene"
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#setwd("D:/GitHub/powerfamily/cross-validation/02_11")
setwd(paste("/home/wang3660/Research/PF/crossvalidation/02_12/", dat, sep=""))

rs = 11
set.seed(rs)
save(rs, file="rs.rda")
load("rs.rda")

load(paste(datapath, dat, ".rda", sep=""))
if(dat == "arcene")
{
  x = tx
  y = ty
  y = c(-1,1)[as.factor(y)]
}

if(dat == "breast")
{
  y = c(-1,1)[as.factor(y)]
}

if(dat == "colon")
{
  x = colon.x
  y = colon.y
  y = c(-1,1)[as.factor(y)]
}

if(dat == "leuk")
{
  y = c(-1,1)[as.factor(y)]
}

if(dat == "prostate")
{
  x = prostate.x
  y = prostate.y
  y = c(-1,1)[as.factor(y)]
}

x = as.matrix(x)
y = c(-1,1)[as.factor(y)]

lambda2seq = c(0,1e-4,1e-3,1e-2,1e-1,1,5,10)
qvseq = c(0.5, 1, 2, 5, 100)

set.seed(rs)
index = sample(1:nrow(x),as.integer(nrow(x)/3),replace=F) 
test_x = x[index,]
test_y = y[index]
train_x = x[-index,]
train_y = y[-index]


if.std = 1
if(if.std == 1)
{
  train_m = apply(train_x, 2, mean)
  std_train_x = t(apply(train_x, 1, function(x) x - train_m))  
  train_sd = apply(std_train_x, 2, function(x) sqrt(x %*% x / length(x)))
  train_sd[train_sd==0] = 1
  std_train_x = t(apply(std_train_x, 1, function(x) x / train_sd))  
  
  std_test_x = t(apply(test_x, 1, function(x) x - train_m))  
  std_test_x = t(apply(std_test_x, 1, function(x) x / train_sd)) 
  
  train_x = std_train_x
  test_x = std_test_x
}

nrep = length(qvseq)
ans = matrix(0, length(lambda2seq), nrep)
l1 = matrix(0, length(lambda2seq), nrep)
tim = matrix(0, length(lambda2seq), nrep)
nzo = matrix(0, length(lambda2seq), nrep)

j = 1
#qv=100
for(qv in qvseq)
{
  i = 1
  print(i)
  for(lambda2 in lambda2seq){
    set.seed(rs)
    fittime = system.time(cv<-cv.GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                      lambda2=lambda2, method="power",
                                      standardize = FALSE,
                                      pred.loss="misclass", nfolds=5))[3]
    cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
    coef.cvs = coef(cv, s="lambda.1se")[-1,]
    nzov = length(coef.cvs[coef.cvs != 0])
    ans[i,j] = cvm
    l1[i,j] = cv$lambda.1se
    tim[i,j] = fittime
    nzo[i,j] = nzov
    i = i + 1
  }
  j = j + 1
}


rownames(ans)=lambda2seq
colnames(ans)=c(qvseq)
rownames(l1)=lambda2seq
colnames(l1)=c(qvseq)
rownames(nzo)=lambda2seq
colnames(nzo)=c(qvseq)
rownames(tim)=lambda2seq
colnames(tim)=c(qvseq)

save(ans, file=paste("ans_", dat, ".rda", sep=""))
save(l1, file=paste("l1_", dat, ".rda", sep=""))
save(nzo, file=paste("nzo_", dat, ".rda", sep=""))
save(tim, file=paste("tim_", dat, ".rda", sep=""))

ans
l1
nzo
tim



bestl2 = apply(ans, 2, function(x) length(lambda2seq) - which.min(rev(x)) + 1 )
cvl2 = lambda2seq[bestl2]
cvl1 = l1[cbind(bestl2, 1:5)]
cvnzo = nzo[cbind(bestl2, 1:5)]
cvans = ans[cbind(bestl2, 1:5)]
cvtim = tim[cbind(bestl2, 1:5)]

bestl2
cvl2
cvl1
cvnzo
cvans
cvtim


save(bestl2, file=paste("bestl2_", dat, ".rda", sep=""))
save(cvl2, file=paste("cvl2_", dat, ".rda", sep=""))
save(cvl1, file=paste("cvl1_", dat, ".rda", sep=""))
save(cvnzo, file=paste("cvnzo_", dat, ".rda", sep=""))
save(cvans, file=paste("cvans_", dat, ".rda", sep=""))
save(cvtim, file=paste("cvtim_", dat, ".rda", sep=""))


