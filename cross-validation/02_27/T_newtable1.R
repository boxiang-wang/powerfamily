rm(list=ls(all=T))
# source("D:/GitHub/powerfamily/U_tool.R")
# datapath = "D:/GitHub/powerfamily/data/"
# dyn.load("D:/GitHub/powerfamily/O_hsvmlassoNET.dll")

source("/home/wang3660/Research/PF/U_tool_server.R")
datapath = "/home/wang3660/Research/PF/data/"
dyn.load("/home/wang3660/Research/PF/crossvalidation/O_hsvmlassoNET.so")

rs = 2  ## specify the random seed
dat = "colon"
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#setwd("D:/GitHub/powerfamily/cross-validation/02_11")
setwd(paste("/home/wang3660/Research/PF/crossvalidation/02_27/", dat, sep=""))


require(gcdnet, Sys.getenv("R_LIBS_USER"))


###############################################
## Step 1 Split the data

#save(rs, file="rs.rda")
#load("rs.rda")
set.seed(rs)


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

lambda2seq = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10)
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

###############################################
## Step 2 Tuning Parameters on Training Sets
##        For the Power Family

nrep = length(qvseq)
ans = matrix(0, length(lambda2seq), nrep)
l1 = matrix(0, length(lambda2seq), nrep)
tim = matrix(0, length(lambda2seq), nrep)
nzo = matrix(0, length(lambda2seq), nrep)

start1 = Sys.time()
j = 1
#qv=100
for(qv in qvseq)
{
  i = 1
  print(i)
  for(lambda2 in lambda2seq){
    set.seed(rs)
    fittime = as.numeric(system.time(cv<-cv.GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                      lambda2=lambda2, method="power",
                                      standardize = FALSE,
                                      pred.loss="misclass", nfolds=5))[3])
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
end1 = Sys.time()

# timepower is the total time of searching the best triple
(timepower = difftime(end1, start1, units="sec"))

# These values are the best values for each fixed q
bestl2_power = apply(ans, 2, function(x) length(lambda2seq) - which.min(rev(x)) + 1 )
cvl2 = lambda2seq[bestl2_power]
cvl1 = l1[cbind(bestl2_power, 1:5)]
cvnzo = nzo[cbind(bestl2_power, 1:5)]
cvans = ans[cbind(bestl2_power, 1:5)]
cvtim = tim[cbind(bestl2_power, 1:5)]

# Now, we select the best q by cv prediction error
tt = cbind(qvseq, cvl2, cvl1, cvans, cvnzo, cvtim, 1:5)
tb1 = tt[order(tt[,5]),]
(ind = as.numeric(tb1[which.min(tb1[,4]), 7]))


###############################################
## Step 3 Tuning Parameters on Training Sets
##        For the logistic regression

ans1 = rep(NA, (length(lambda2seq)))
l11  = rep(NA, (length(lambda2seq))) 
tim1 = rep(NA, (length(lambda2seq)))
nzo1 = rep(NA, (length(lambda2seq)))

# train a model
start2 = Sys.time()
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
end2 = Sys.time()

# timepower is the total time of searching the best triple
timelogis = difftime(end2, start2, units="sec")

bestl2_logis = length(lambda2seq) - which.min(rev(ans1)) + 1
cvl2_logis = lambda2seq[bestl2_logis]
cvl1_logis = l11[bestl2_logis]


###############################################
## Step 4 Test Error


## For the Power Family

ans_q = rep(NA, (length(qvseq)))
tim_q = rep(NA, (length(qvseq)))
nzo_q = rep(NA, (length(qvseq)))

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


## For logistic regression

tims = as.numeric(system.time(fitlog<-gcdnet(train_x, train_y, eps=1e-8, 
                                  lambda2=cvl2_logis, method="logit", 
                                  standardize=F))[3])
fit.coef = as.vector(coef(fitlog, s=cvl1_logis))[-1]
nzov= length(fit.coef[fit.coef != 0]) 
pre = predict(fitlog, newx=test_x, s=cvl1_logis, type="class" )
(ans_logis = mean(test_y != pre) )
(tim_logis = tims)
(nzo_logis = nzov)


## For SVM

i = length(qvseq) + 2
tim1 = as.numeric(system.time(fit1<-tryCatch(DrSVM_Fix2(train_x, train_y, 
                                    eps=1e-8,
                                    lambda2=cvl2[length(qvseq)]*nrow(train_x),
                                    scale=F),
                                    error = function(e) NA))[3])
if(!is.na(fit1))
{
  fit.coef1 = M_coef.powerfamily(newobj(fit1), s=cvl1[length(qvseq)]*nrow(train_x))
  f = fit.coef1[-1]
  f1 = fit.coef1[1]
  nzov1= length(f[f != 0])
  pre1 = sign(test_x %*% f + 
                matrix(rep(f1, nrow(test_x)), nrow=nrow(test_x), byrow=T))
  (ans_SVM = mean(test_y != pre1) )
  (nzo_SVM = nzov1)
  (tim_SVM = tim1)
}else{
  ans_SVM = NA
  nzo_SVM = NA
  tim_SVM = NA
}


(prop.q = qvseq[ind])
(ans.vec = c(ans_q[ind], ans_q, ans_SVM, ans_logis))
(tim.vec = c(tim_q[ind], tim_q, tim_SVM, tim_logis))
save(prop.q, ans.vec, tim.vec, file=paste("res_", dat, rs, ".rda", sep=""))
