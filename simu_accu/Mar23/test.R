tt = 2
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

qv100l2.list[tt] = bestl2
qv100l1.list[tt] = bestl1
mDr = DrSVM_Fix2(train$x, train$y, lambda2=qv100l2.list[tt]*nrow(train$x),
                 eps=1e-8,scale=F, quiet=F)
fit.coef1 = M_coef.powerfamily(newobj(mDr), s=qv100l1.list[tt]*nrow(train$x))
f = fit.coef1[-1]
f1 = fit.coef1[1]
pre1 = sign(test$x %*% f + matrix(rep(f1, nrow(test$x)), nrow=nrow(test$x), byrow=T))
error_SVM[tt] = mean(test$y != pre1)
C_SVM[tt] = sum(f[1:p0] != 0)
I_SVM[tt] = sum(f[-c(1:p0)] != 0)

mDr = DrSVM_Fix2(train$x, train$y, lambda2=10*nrow(train$x),
             eps=1e-8,scale=F, quiet=F)  