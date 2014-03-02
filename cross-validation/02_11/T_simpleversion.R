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


#setwd("D:/GitHub/powerfamily/cross-validation/02_11")
setwd(paste("/home/wang3660/Research/PF/crossvalidation/02_12/", dat, sep=""))


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


qv=5
lambda2=5
j = 1
for(qv in qvseq)
{
  i = 1
  print(i)
  for(lambda2 in lambda2seq){
    set.seed(1234)
    fittime = system.time(cv<-cv.GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                      lambda2=lambda2, method="power",
                                      standardize=F,
                                      pred.loss="misclass", nfolds=5))[3]
    cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
    coef.cvs = coef(cv, s="lambda.1se")[-1,]
    nzov = length(coef.cvs[coef.cvs != 0]) 
    # pre = predict(cv$GCDpower.fit, newx = test_x, s = cv$lambda.1se, type = "class" )
    # error = (test_y != pre) 
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

# load(paste("ans_", dat, ".rda", sep="")); load(paste("l1_", dat, ".rda", sep=""))
# load(paste("nzo_", dat, ".rda", sep="")); load(paste("tim_", dat, ".rda", sep=""))

bestl2 = apply(ans, 2, function(x) length(lambda2seq) - which.min(rev(x)) + 1 )
cvl2 = lambda2seq[bestl2]
cvl1 = l1[cbind(bestl2, 1:5)]
cvnzo = nzo[cbind(bestl2, 1:5)]
cvans = ans[cbind(bestl2, 1:5)]
 

ans_q = rep(NA, (length(qvseq) + 2))
tim_q = rep(NA, (length(qvseq) + 2))
nzo_q = rep(NA, (length(qvseq) + 2))

i = length(qvseq) 
for(i in 1:length(qvseq))
{
  print(i)
  tims = system.time(fit<-GCDpower(train_x, train_y, eps=1e-8, qv=qvseq[i],
                                    lambda2=cvl2[i], method="power", 
                                   standardize=F))[3]
  fit.coef = as.vector(coef(fit, s=cvl1[i]))
  nzov= length(fit.coef[fit.coef != 0]) - 1
  pre = predict(fit, newx=test_x, s=cvl1[i], type="class" )
  ans_q[i] = mean(test_y != pre) 
  nzo_q[i] = nzov
  tim_q[i] = tims
}

i = length(qvseq) + 1
print(i)
tim = system.time(fith<-GCDpower(train_x, train_y, eps=1e-8, delta=0.01,
                                lambda2=cvl2[(i-1)], method="hhsvm"))[3]
fith.coef = as.vector(coef(fith, s=cvl1[(i-1)]))
nzov= length(fith.coef[fith.coef != 0]) - 1
pre = predict(fith, newx=test_x, s=cvl1[(i-1)], type="class" )
ans_q[i] = mean(test_y != pre) 
nzo_q[i] = nzov
tim_q[i] = tim

fit10 <- GCDpower(train_x, train_y, eps=1e-8, qv=10, lambda2=1, method="power", standardize=F)
fit100 <- GCDpower(train_x, train_y, eps=1e-8, qv=100, lambda2=1, method="power", standardize=F)
fit101 <- GCDpower(train_x, train_y, eps=1e-8, qv=101, lambda2=1, method="power", standardize=F)
fit001 <- GCDpower(train_x, train_y, eps=1e-8, delta=0.01, lambda2=1, method="hhsvm", standardize=F)

fitd <- DrSVM_Fix2(train_x, train_y, eps=1e-8, lambda2=1*nrow(train_x), scale=F)

xlim=c(0,8); ylim=c(-0.18, 0.15)
plot(fit100
     , xlim=xlim,  ylim=ylim
     )
plot(fit001
     , xlim=xlim,  ylim=ylim
     )








cc1 = M_coef.powerfamily(newobj(fitd), s = (fit100$lambda*nrow(train_x)))[-1,]
dimnames(cc1)=NULL
max(abs(fit100$beta - cc1))
max(abs(fit001$beta - cc1))
max(abs(fit100$beta-fit001$beta))


cc = M_coef.powerfamily(newobj(fitd), s = fit100$lambda*nrow(train_x))
s = apply(abs(cc[-1,]), 2, sum)
mm = cbind(t(cc[-1,]))
par(mar=c(5.1,5.1,4.1,2.1))
matplot(s, mm, type="l", cex.lab = 1.5, 
        lty=1, col = gray.colors(12, start = 0.05,
                                 end = 0.7, gamma = 2.2),
        
        xlim=c(0,8), ylim=c(-0.18, 0.15),
        xlab="L1 Norm", 
        ylab="Coefficients")




ans_q
nzo_q
tim_q

i = length(qvseq) + 2
tim1 = system.time(fit1<-DrSVM_Fix2(train_x, train_y, 
                                  eps=1e-8, lambda2=cvl2[length(qvseq)]*nrow(train_x),
                                  scale=F))[3]
fit.coef1 = M_coef.powerfamily(newobj(fit1), s=cvl1[length(qvseq)]*nrow(train_x))
nzov1= length(fit.coef1[fit.coef1 != 0]) - 1

pre1 = sign(test_x %*% fit.coef1[-1,] + 
       matrix(rep(fit.coef1[1,], nrow(test_x)), nrow=nrow(test_x), byrow=T))
ans_q[i] = mean(test_y != pre1) 
nzo_q[i] = nzov1
tim_q[i] = tim1

#ans_q
#nzo_q
#tim_q
#cvl1
#cvl2

#res = cbind(ans_q, nzo_q, tim_q, cvl1, cvl2)

#save(res, file=paste("res_", dat, ".rda", sep=""))




