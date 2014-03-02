# source("D:/GitHub/powerfamily/U_tool.R")
# datapath = "D:/GitHub/powerfamily/data/"
source("/home/wang3660/Research/PF/U_tool_server.R")
datapath = "/home/wang3660/Research/PF/data/"

# dyn.load("D:/GitHub/powerfamily/O_hsvmlassoNET.dll")
dyn.load("/home/wang3660/Research/PF/crossvalidation/O_hsvmlassoNET.so")

dat = "colon"
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
setwd(paste("/home/wang3660/Research/PF/crossvalidation/02_12/", dat, sep=""))
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

load(paste("ans_", dat, ".rda", sep=""))
load(paste("l1_", dat, ".rda", sep=""))
load(paste("nzo_", dat, ".rda", sep=""))
load(paste("tim_", dat, ".rda", sep=""))


load(paste("bestl2_", dat, ".rda", sep=""))
load(paste("cvl2_", dat, ".rda", sep=""))
load(paste("cvl1_", dat, ".rda", sep=""))
load(paste("cvnzo_", dat, ".rda", sep=""))
load(paste("cvans_", dat, ".rda", sep=""))
load(paste("cvtim_", dat, ".rda", sep=""))

load(paste("ans_q_2", dat, ".rda", sep="")); 
load(paste("nzo_q_2", dat, ".rda", sep="")); 
load(paste("tim_q_2", dat, ".rda", sep=""))

## Create the personal library if it doesn't exist. Ignore a warning 
# if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
## Install one package.
install.packages("gcdnet", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
require(gcdnet, Sys.getenv("R_LIBS_USER"))





##############################
# logistic
ans1 = rep(NA, (length(lambda2seq)))
l11  = rep(NA, (length(lambda2seq))) 
tim1 = rep(NA, (length(lambda2seq)))
nzo1 = rep(NA, (length(lambda2seq)))

# train a model
i = 1
for(lambda2 in lambda2seq){
  set.seed(rs)
  fittime = system.time(cv<-cv.gcdnet(train_x, train_y, eps=1e-8, 
                                        lambda2=lambda2, method="logit",
                                        standardize = FALSE,
                                        pred.loss="misclass", nfolds=5))[3]
  cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
  coef.cvs = coef(cv, s="lambda.1se")[-1,]
  nzov = length(coef.cvs[coef.cvs != 0])
  ans1[i] = cvm
  l11[i] = cv$lambda.1se
  tim1[i] = fittime
  nzo1[i] = nzov
  i = i + 1
}


ans1
nzo1
l11
tim1

bestl2 = length(lambda2seq) - which.min(rev(ans1)) + 1

# test a model
tims = system.time(fitlog<-gcdnet(train_x, train_y, eps=1e-8, 
                                 lambda2=lambda2seq[bestl2], method="logit", 
                                 standardize=F))[3]
fit.coef = as.vector(coef(fitlog, s=l11[bestl2]))[-1]
nzov= length(fit.coef[fit.coef != 0]) 
pre = predict(fitlog, newx=test_x, s=l11[bestl2], type="class" )
(ans_logis = mean(test_y != pre) )
(tim_logis = tims)
(nzo_logis = nzov)

#########################
# sqsvm

ans2 = rep(NA, (length(lambda2seq)))
l12  = rep(NA, (length(lambda2seq))) 
tim2 = rep(NA, (length(lambda2seq)))
nzo2 = rep(NA, (length(lambda2seq)))

# train a model
i = 1
for(lambda2 in lambda2seq){
  set.seed(rs)
  fittime = system.time(cv<-cv.gcdnet(train_x, train_y, eps=1e-8, 
                                      lambda2=lambda2, method="sqsvm",
                                      standardize = FALSE,
                                      pred.loss="misclass", nfolds=5))[3]
  cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
  coef.cvs = coef(cv, s="lambda.1se")[-1,]
  nzov = length(coef.cvs[coef.cvs != 0])
  ans2[i] = cvm
  l12[i] = cv$lambda.1se
  tim2[i] = fittime
  nzo2[i] = nzov
  i = i + 1
}


ans2
nzo2
l12
tim2

bestl2 = length(lambda2seq) - which.min(rev(ans1)) + 1

# test a model
tims = system.time(fitlog<-gcdnet(train_x, train_y, eps=1e-8, 
                                  lambda2=lambda2seq[bestl2], method="sqsvm", 
                                  standardize=F))[3]
fit.coef = as.vector(coef(fitlog, s=l12[bestl2]))[-1]
nzov= length(fit.coef[fit.coef != 0]) 
pre = predict(fitlog, newx=test_x, s=l12[bestl2], type="class" )
(ans_sqsvm = mean(test_y != pre) )
(tim_sqsvm = tims)
(nzo_sqsvm = nzov)


#########################
# hhsvm (delta = 2)

ans3 = rep(NA, (length(lambda2seq)))
l13  = rep(NA, (length(lambda2seq))) 
tim3 = rep(NA, (length(lambda2seq)))
nzo3 = rep(NA, (length(lambda2seq)))

# train a model
i = 1
for(lambda2 in lambda2seq){
  set.seed(rs)
  fittime = system.time(cv<-cv.gcdnet(train_x, train_y, eps=1e-8, delta=2,
                                      lambda2=lambda2, method="hhsvm",
                                      standardize = FALSE,
                                      pred.loss="misclass", nfolds=5))[3]
  cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
  coef.cvs = coef(cv, s="lambda.1se")[-1,]
  nzov = length(coef.cvs[coef.cvs != 0])
  ans3[i] = cvm
  l13[i] = cv$lambda.1se
  tim3[i] = fittime
  nzo3[i] = nzov
  i = i + 1
}


ans3
nzo3
l13
tim3

bestl2 = length(lambda2seq) - which.min(rev(ans1)) + 1

# test a model
tims = system.time(fitlog<-gcdnet(train_x, train_y, eps=1e-8, delta=2,
                                  lambda2=lambda2seq[bestl2], method="hhsvm", 
                                  standardize=F))[3]
fit.coef = as.vector(coef(fitlog, s=l13[bestl2]))[-1]
nzov= length(fit.coef[fit.coef != 0]) 
pre = predict(fitlog, newx=test_x, s=l13[bestl2], type="class" )
(ans_hhsvm = mean(test_y != pre) )
(tim_hhsvm = tims)
(nzo_hhsvm = nzov)


save(ans_sqsvm, tim_sqsvm, nzo_sqsvm, ans_logis, tim_logis, nzo_logis, ans_hhsvm, 
     tim_hhsvm, nzo_hhsvm, file=paste("others_", dat, ".rda", sep=""))

  load(paste("ans_q_2", dat, ".rda", sep=""))
  load(paste("nzo_q_2", dat, ".rda", sep=""))
  load(paste("tim_q_2", dat, ".rda", sep=""))
  
  
  load(paste("ans_ql1_", dat, ".rda", sep=""))
  load(paste("nzo_ql1_", dat, ".rda", sep=""))
  load(paste("tim_ql1_", dat, ".rda", sep=""))


tt = cbind(qvseq, cvl2, cvl1, cvans, cvnzo, cvtim, 1:5)
tb1 = tt[order(tt[,5]),]
(ind = as.numeric(tb1[which.min(tb1[,4]), 7]))

(aa = c(qvseq[ind], ans_q[ind], nzo_q[ind], tim_q[ind], ans_q[7], nzo_q[7], tim_q[7],
  ans_sqsvm, nzo_sqsvm, tim_sqsvm, ans_logis, nzo_logis, tim_logis, 
  ans_hhsvm, nzo_hhsvm, tim_hhsvm, ans_ql1[2], nzo_ql1[2], tim_ql1[2]))
save(aa, file=paste("res_", dat, ".rda", sep=""))


aa2 = c(qvseq[ind],cvl2[ind], cvl1[ind], cvans[ind], cvtim[ind])
save(aa2, file=paste("train_", dat, ".rda", sep=""))

  nzo
  ans
  bestl2
  cvans
  cvnzo
  cvl2
  cvl1
  ans_q
  nzo_q
  tim_q
  ans_ql1
  nzo_ql1
  tim_ql1
  ans_sqsvm
  tim_sqsvm
  nzo_sqsvm
  
  ans_logis
  tim_logis
  nzo_logis

  ans_hhsvm
  tim_hhsvm
  nzo_hhsvm



