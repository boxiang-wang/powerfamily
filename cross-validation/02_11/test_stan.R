 source("D:/GitHub/powerfamily/U_tool.R")
 datapath = "D:/GitHub/powerfamily/data/"
#source("/home/wang3660/Research/PF/U_tool_server.R")
#datapath = "/home/wang3660/Research/PF/data/"

# dyn.load("D:/GitHub/powerfamily/O_hsvmlassoNET.dll")
#dyn.load("/home/wang3660/Research/PF/crossvalidation/O_hsvmlassoNET.so")

dat = "breast"
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


setwd("D:/GitHub/powerfamily/cross-validation/02_11")
# setwd(paste("/home/wang3660/Research/PF/crossvalidation/02_12/", dat, sep=""))


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

set.seed(1234)
index = sample(1:nrow(x),as.integer(nrow(x)/5),replace=F) 
test_x = x[index,]
test_y = y[index]
train_x = x[-index,]
train_y = y[-index]

 qv = 1
 l2 = 5
 
 (tims1 = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                   lambda2=l2, method="power", 
                                   standardize=F))[3])
 fit.coef = as.vector(coef(fit, s=0.01))
 (nzov1= length(fit.coef[fit.coef != 0]) - 1)
 (pre = as.vector(predict(fit, newx=test_x, s=0.01, type="class" )))
 (e1 = mean(test_y != pre) )
 
 (tims2 = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                    lambda2=l2, method="power", 
                                    standardize=T))[3])
 fit.coef = as.vector(coef(fit, s=0.01))
 (nzov2= length(fit.coef[fit.coef != 0]) - 1)
 (pre = as.vector(predict(fit, newx=test_x, s=0.01, type="class" )))
 (e2 = mean(test_y != pre) )


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


 
 (tims3 = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                    lambda2=l2, method="power", 
                                    standardize=F))[3])
 fit.coef = as.vector(coef(fit, s=0.01))
 (nzov3= length(fit.coef[fit.coef != 0]) - 1)
 (pre = as.vector(predict(fit, newx=test_x, s=0.01, type="class" )))
 (e3 = mean(test_y != pre) )
 
 (tims4 = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                    lambda2=l2, method="power", 
                                    standardize=T))[3])
 fit.coef = as.vector(coef(fit, s=0.01))
 (nzov4= length(fit.coef[fit.coef != 0]) - 1)
 (pre = as.vector(predict(fit, newx=test_x, s=0.01, type="class" )))
 (e4 = mean(test_y != pre) )
 