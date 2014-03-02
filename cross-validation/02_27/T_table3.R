# source("D:/GitHub/powerfamily/U_tool.R")
# datapath = "D:/GitHub/powerfamily/data/"
source("/home/wang3660/Research/PF/U_tool_server.R")
datapath = "/home/wang3660/Research/PF/data/"

# dyn.load("D:/GitHub/powerfamily/O_hsvmlassoNET.dll")
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


load(paste("ans_q_", dat, ".rda", sep=""))
load(paste("nzo_q_", dat, ".rda", sep=""))
load(paste("tim_q_", dat, ".rda", sep=""))

i = length(qvseq) + 2
tim1 = system.time(fit1<-DrSVM_Fix2(train_x, train_y, 
                                    eps=1e-8, 
                                    lambda2=cvl2[length(qvseq)]*nrow(train_x),
                                    scale=F))[3]
fit.coef1 = M_coef.powerfamily(newobj(fit1), s=cvl1[length(qvseq)]*nrow(train_x))
f = fit.coef1[-1]
f1 = fit.coef1[1]

nzov1= length(f[f != 0])

pre1 = sign(test_x %*% f + 
              matrix(rep(f1, nrow(test_x)), nrow=nrow(test_x), byrow=T))
ans_q[i] = mean(test_y != pre1) 
nzo_q[i] = nzov1
tim_q[i] = tim1

ans_q
nzo_q
tim_q


save(ans_q, file=paste("ans_q_2", dat, ".rda", sep=""))
save(nzo_q, file=paste("nzo_q_2", dat, ".rda", sep=""))
save(tim_q, file=paste("tim_q_2", dat, ".rda", sep=""))


if(1 == 2){
  
dat = "prostate"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_12/", dat, sep=""))
load(paste("ans_q_2", dat, ".rda", sep=""))
load(paste("nzo_q_2", dat, ".rda", sep=""))
load(paste("tim_q_2", dat, ".rda", sep=""))

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
}