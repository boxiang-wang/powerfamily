rm(list=ls(all=TRUE))

args=(commandArgs(TRUE))
if(length(args)==0){
 print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
   eval(parse(text=args[[i]]))
  }
}
 
if(!exists("nms")) nms=12
print(nms) 
 
if(.Platform$OS.type == "unix")
{
  file.loc = paste("/home/wang3660/Research/PF/simu_accu/Mar23/Mar25res/", nms, '/', sep="")
  setwd(file.loc)
  source("/home/wang3660/Research/PF/U_tool_server.R")
  require(gcdnet, Sys.getenv("R_LIBS_USER"))
}


if(.Platform$OS.type == "windows")
{
  file.loc = paste("D:/GitHub/powerfamily_output/simu_accu/Mar25res/", nms, '/', sep="")
  setwd(file.loc)
  source("D:/GitHub/powerfamily/U_tool.R")
  require(gcdnet)
}


# s_tr is the size of training set.
# s_vl is the size of validation set, which is used to find tuning parameters.
# s_ts is the size of test set, which is used to compare methods.

if(!exists("s_tr")) s_tr=100
if(!exists("s_vl")) s_vl=100
if(!exists("s_ts")) s_ts=30000
print(s_tr); print(s_vl); print(s_ts); 
 
# p0 is the number of relevant variables.
if(!exists("pp")) pp=300
if(!exists("p0")) p0=5
print(pp); print(p0);
 
# mu is the mean shift
if(!exists("mu")) mu=0.6
print(mu);

# rho is the correlation
if(!exists("rho")) rho=0.8
print(rho);
 
# blanc is the size ratio of plus and minus groups
if(!exists("blanc")) blanc=0.5
print(blanc);

# senario of covariance matrix
if(!exists("sina")) sina=2
file.loc = paste(file.loc, 'sina', sina, sep="")

# whether include DrSVM
if(!exists("inc.Dr")) inc.Dr=T
if(inc.Dr == T) file.loc = paste(file.loc, "Dr", sep="")
setwd(file.loc)

# simulation times
if(!exists("simu")) simu=5
print(simu);

# start and end points of the simulation 
 
if(!exists('st')) st = 1
if(!exists('en')) en = simu
print(st); print(en) 
 
print(file.loc) 
# Bign is the total number of generated observations.
Bign = s_tr + s_vl + s_ts
Bigx = matrix(NA, Bign, pp)
Bigy = rep(NA, Bign)

# The means of two classes.
mup = c(rep(mu, p0), rep(0, pp - p0))
mum = c(rep(-mu, p0), rep(0, pp - p0))

# The var-cov matrix of two classes.
if(sina == 1){
  sigma = diag(pp)
}
if(sina == 2){
  sigma = diag(pp)
  for(i in 1:p0) for(j in 1:p0) 
    sigma[i,j] = ifelse(i == j, 1, rho)
  
}
if(sina == 3){
  sigma = diag(pp)
  for(i in 1:p0) for(j in 1:p0) 
    sigma[i,j] = rho ^ abs(i-j)
}

eo = eigen(sigma, symmetric=TRUE)
sigma.sqrt = eo$vec %*% (diag(sqrt(eo$val)) %*% t(eo$vec))

# The index of two classes.
plusgroup = c(1:(s_tr*blanc), s_tr + 1:(s_vl*blanc), s_tr + s_vl + 1:(s_ts*blanc))
minusgroup = (1:Bign)[-plusgroup]

l2.list = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10)
qv.list = c(0.5, 1, 2, 5, 100)
l2.length = length(l2.list); qv.length = length(qv.list)


if(st == 1){
  # For elastic net and l1, they record the prediction error and C, I.
  qverror_elas = qverror_l1 = matrix(NA, simu, qv.length)
  C_elas = C_l1 = matrix(NA, simu, qv.length)
  I_elas = I_l1 = matrix(NA, simu, qv.length)
  
  # They record the best q tuned by cross validation
  qverror_cv_elas = C_cv_elas = I_cv_elas = rep(NA, simu)
  qverror_cv_l1 = C_cv_l1 = I_cv_l1 = rep(NA, simu)
  
  # They record the error for SVM and the best lambda for q=100
  if(inc.Dr == T){
    error_SVM = C_SVM = I_SVM = rep(NA, simu)
    qv100l2.list = rep(NA, simu)
    qv100l1.list = rep(NA, simu)
  }
} else {
    load("allres.rda")
    if(inc.Dr == T){
      load("allSVM.rda")
    }
}
 

set.seed(1234)
seeds = sample(1:10000, size=10000)
for(tt in st:en){
  print(tt)
  write.csv(tt, file=paste(tt, ".csv", sep=""))
  set.seed(seeds[tt])
  for(i in plusgroup){
    Bigx[i,] = mup + sigma.sqrt %*% rnorm(pp)
    Bigy[i] = 1
  }
  for( i in minusgroup){
    Bigx[i,] = mum + sigma.sqrt %*% rnorm(pp)
    Bigy[i] = -1 
  }
  
  train = list(x = Bigx[1:s_tr, ], y = Bigy[1:s_tr])
  validation = list(x = Bigx[s_tr + 1:s_vl, ], y = Bigy[s_tr + 1:s_vl])
  test = list(x = Bigx[s_tr + s_vl + 1:s_ts, ], y = Bigy[s_tr + s_vl + 1:s_ts])
  
  # standardize the input
  # both the validation and the testing set are standardized by the training set
  train_m = apply(train$x, 2, mean)
  std_train_x = t(apply(train$x, 1, function(x) x - train_m))  
  train_sd = apply(std_train_x, 2, function(x) sqrt(x %*% x / length(x)))
  train_sd[train_sd==0] = 1
  std_train_x = t(apply(std_train_x, 1, function(x) x / train_sd))  
  
  std_test_x = t(apply(test$x, 1, function(x) x - train_m))  
  std_test_x = t(apply(std_test_x, 1, function(x) x / train_sd)) 
  
  std_vali_x = t(apply(validation$x, 1, function(x) x - train_m))  
  std_vali_x = t(apply(std_vali_x, 1, function(x) x / train_sd)) 
  
  train$x = std_train_x
  test$x = std_test_x
  validation$x = std_vali_x
  
  # qvcverror records the minimum error for each q
  qvcverror = rep(NA, qv.length)
  qvcvnzo = rep(NA, qv.length)
  
  
  for(i in 1:qv.length){
    qv = qv.list[i]
    l2error = rep(NA, l2.length)
    l2nzo = rep(NA, l2.length)
    l2bestl1 = rep(NA, l2.length)
    modellist = as.list(1:l2.length)
    for (j in 1:l2.length){
      # fit models on the training set
      m = GCDpower(train$x, train$y, lambda2=l2.list[j], qv=qv.list[i], 
                   method="power", eps=1e-8, strong=T, standardize=T) 
      
      # select tuning parameters on the validation set
      pre = predict(m, type="class", newx=validation$x)
      vali.error = apply(pre, 2, function(cols) mean(cols != validation$y))
      l2error[j] = min(vali.error) 
      l2bestl1[j] = m$lambda[which.min(vali.error)]
      l2nzo[j] = length(coef(m, type="nonzero", s=l2bestl1[j]))
      modellist[[j]] = m
    }
    
    qvcverror[i] = min(l2error)
    qvcvnzo[i] = l2nzo[which.min(l2error)]
    
    # find the best tuning paramter
    bestidl2 = tail(which(l2error == min(l2error)) , 1)
    bestl2 = l2.list[bestidl2]
    bestl1 = l2bestl1[bestidl2]
    
    # for elastic net penalty:
    # compute the accuracy rate on the test set
    pre = predict(modellist[[bestidl2]], type="class", newx=test$x, s=bestl1)
    qverror_elas[tt, i] = mean(pre != test$y)
    m.coef = coef(modellist[[bestidl2]], type="nonzero", s=bestl1)
    C_elas[tt, i] = sum(c(1:p0) %in% m.coef)
    I_elas[tt, i] = length(m.coef) - C_elas[tt, i]
       
    # for l1 penalty
    pre1 = predict(modellist[[1]], type="class", newx=test$x, s=l2bestl1[1])
    qverror_l1[tt, i] = mean(pre1 != test$y)
    m.coef1 = coef(modellist[[1]], type="nonzero", s=l2bestl1[1])
    C_l1[tt, i] = sum(c(1:p0) %in% m.coef1)
    I_l1[tt, i] = length(m.coef1) - C_l1[tt, i]
    
  }
  
  if(inc.Dr == T){
    
    qv100l2.list[tt] = bestl2
    qv100l1.list[tt] = bestl1
    
    l2Dr = max(qv100l2.list[tt], 1e-6)
    mDr <- tryCatch(DrSVM_Fix2(train$x, train$y, lambda2=l2Dr*nrow(train$x),
                               eps=1e-8,scale=F), 
                    error = function(e) NA)
    if(!is.na(mDr)){
      fit.coef1 = M_coef.powerfamily(newobj(mDr), s=qv100l1.list[tt]*nrow(train$x))
      f = fit.coef1[-1]
      f1 = fit.coef1[1]
      pre1 = sign(test$x %*% f + matrix(rep(f1, nrow(test$x)), nrow=nrow(test$x), byrow=T))
      error_SVM[tt] = mean(test$y != pre1)
      C_SVM[tt] = sum(f[1:p0] != 0)
      I_SVM[tt] = sum(f[-c(1:p0)] != 0)
    } else {
      error_SVM[tt] = NA
      C_SVM[tt] = NA
      I_SVM[tt] = NA
    }
    
    
  }
  
  # find the best q by cv
  toget = cbind(qvcverror, qvcvnzo, 1:qv.length)
  toget1 = toget[order(toget[,2]),]
  qvbycv_ind = as.numeric(toget1[which.min(toget1[,1]), 3])
  
  qverror_cv_elas[tt] = qverror_elas[tt, qvbycv_ind]
  C_cv_elas[tt] = C_elas[tt, qvbycv_ind]
  I_cv_elas[tt] = I_elas[tt, qvbycv_ind]
  
  qverror_cv_l1[tt] = qverror_l1[tt, qvbycv_ind]
  C_cv_l1[tt] = C_l1[tt, qvbycv_ind]
  I_cv_l1[tt] = I_l1[tt, qvbycv_ind]
  
}

qverror_l1 = cbind(qverror_l1, qverror_cv_l1)
qverror_elas = cbind(qverror_elas, qverror_cv_elas)
C_l1 = cbind(C_l1, C_cv_l1)
C_elas = cbind(C_elas, C_cv_elas)
I_l1 = cbind(I_l1, I_cv_l1)
I_elas = cbind(I_elas, I_cv_elas)

dimnames(qverror_l1) = dimnames(qverror_elas) = NULL
dimnames(C_l1) = dimnames(C_elas) = dimnames(I_l1) = dimnames(I_elas) = NULL


qv.list

apply(qverror_l1, 2, median)
apply(qverror_elas, 2, median)

apply(C_l1, 2, median)
apply(C_elas, 2, median)

apply(I_l1, 2, median)
apply(I_elas, 2, median)

save(qverror_l1, qverror_elas, file="qverror.rda")
save(C_l1, C_elas, I_l1, I_elas, file="CandI.rda")

if(inc.Dr == T){
  qverror_elas = cbind(qverror_elas, error_SVM)
  C_elas = cbind(C_elas, C_SVM)
  I_elas = cbind(I_elas, I_SVM)
  
  save(qv100l2.list, qv100l1.list, file="qv100lambda.rda")
}
 
apply(na.omit(qverror_elas), 2, median)
apply(na.omit(C_elas), 2, median)
apply(na.omit(I_elas), 2, median)
 
save(qverror_elas, qverror_l1, C_elas, C_l1, I_elas, I_l1, qverror_cv_elas, C_cv_elas, 
     I_cv_elas, qverror_cv_l1, C_cv_l1, I_cv_l1, file="allres.rda")
if(inc.Dr == T){
  save(error_SVM, C_SVM, I_SVM, qv100l2.list, qv100l1.list, file="allSVM.rda")
}
