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
  id = 1
}

if(dat == "breast")
{
  id = 2
}

if(dat == "colon")
{
  x = colon.x
  y = colon.y
  id = 3
}

if(dat == "leuk")
{
  id = 4
}

if(dat == "prostate")
{
  x = prostate.x
  y = prostate.y
  id = 5
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

load('res.rda')
aa



i = length(qvseq) 
for(i in 1:length(qvseq))
{
  print(i)
  tims = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qvseq[i],
                                   lambda2=cvl2[i], method="power", 
                                   standardize=F))[3]
  fit.coef = as.vector(coef(fit, s=cvl1[i]))[-1]
  nzov= length(fit.coef[fit.coef != 0])
  pre = predict(fit, newx=test_x, s=cvl1[i], type="class" )
  ans_q[i] = mean(test_y != pre) 
  nzo_q[i] = nzov
  tim_q[i] = tims
}





i = length(qvseq) + 1
print(i)
tims = system.time(fith<-GCDpower(train_x, train_y, eps=1e-8, delta=0.01,
                                  lambda2=cvl2[(i-1)], method="hhsvm", standardize=F))[3]
fith.coef = as.vector(coef(fith, s=cvl1[(i-1)]))[-1]
nzov= length(fith.coef[fith.coef != 0])
pre = predict(fith, newx=test_x, s=cvl1[(i-1)], type="class" )
ans_q[i] = mean(test_y != pre) 
nzo_q[i] = nzov
tim_q[i] = tims


save(ans_q, file=paste("ans_q_", dat, ".rda", sep=""))
save(nzo_q, file=paste("nzo_q_", dat, ".rda", sep=""))
save(tim_q, file=paste("tim_q_", dat, ".rda", sep=""))



nzo
ans
bestl2
cvans
cvnzo
cvl2
cvl1
ans_q
nzo_q
