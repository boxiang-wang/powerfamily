# rm(list=ls(all=TRUE))

source("D:\\GitHub\\powerfamily\\U_tool.R")
setwd("D:\\GitHub\\powerfamily\\cross-validation\\02_11")

load("D:\\GitHub\\powerfamily\\data\\colon.rda")
x = colon.x
y = colon.y

x = as.matrix(x)
y = c(-1,1)[as.factor(y)]


args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}



print(qv)
if(!exists('st'))
{
  st = 1
}
if(!exists('nrep'))
{
  nrep = 10
}
if(!exists('en'))
{
  en = nrep
}

# nrep = 1

# set.seed(1234)
# seed = sample((1:nrep) + 2001, nrep)

lambda2seq = c(1e-4,1e-3,1e-2,1e-1,1,5,10)

if(st == 1)
{
  ans = matrix(0, length(lambda2seq), nrep)
  l1 = matrix(0, length(lambda2seq), nrep)
  time = matrix(0, length(lambda2seq), nrep)
  nzo = matrix(0, length(lambda2seq), nrep)
} else
{
  load(paste("q=", qv, "_ans.rda", sep=""))
  load(paste("q=", qv, "_time.rda", sep=""))
  load(paste("q=", qv, "_nzo.rda", sep=""))
}



write.csv(NA, file=paste("q_", qv, "_j_", 0, ".csv", sep=""))
for (j in st:en){
  print(paste("q =", qv, "j =", j))
  write.csv(NA, file=paste("q_", qv, "_j_", j, ".csv", sep=""))
  file.remove(paste("q_", qv, "_j_", (j-1), ".csv", sep=""))
  i = 1
  set.seed(1234)
  index = sample(1:nrow(x),as.integer(nrow(x)/3),replace=F) 
  test_x = x[index,]
  test_y = y[index]
  train_x = x[-index,]
  train_y = y[-index]
  
  train_m = apply(train_x, 2, mean)
  std_train_x = t(apply(train_x, 1, function(x) x - train_m))  
  train_sd = apply(std_train_x, 2, function(x) sqrt(x %*% x))
  std_train_x = t(apply(std_train_x, 1, function(x) x / train_sd))  

  std_test_x = t(apply(test_x, 1, function(x) x - train_m))  
  std_test_x = t(apply(std_test_x, 1, function(x) x / train_sd))  
  
  for(lambda2 in c(1e-4,1e-3,1e-2,1e-1,1,5,10)){
    set.seed(1234)
    tim = system.time(cv<-cv.GCDpower(std_train_x, train_y, eps=1e-8, qv=qv, delta=2,
                                      lambda2=lambda2, method="power",
                                      pred.loss="misclass", nfolds=5))[3]
    l1[i,j] = cv$lambda.1se
    coef.cvs = coef(cv, s="lambda.1se")[-1,]
    nzov= length(coef.cvs[coef.cvs != 0])
    pre = predict(cv$GCDpower.fit, newx = std_test_x, s = cv$lambda.1se, type = "class" )
    error = (test_y != pre) 
    nzo[i,j] = nzov
    time[i,j] = tim
    ans[i,j] = res = mean(error)
    i = i + 1
  }
  save(ans, file=paste("q=", qv, "_ans.rda", sep=""))
  save(l1, file=paste("q=", qv, "_l1.rda", sep=""))
  save(nzo, file=paste("q=", qv, "_nzo.rda", sep=""))
  save(time, file=paste("q=", qv, "_time.rda", sep=""))
}
(ans.avg = apply(ans,1,mean)*100)
(l1.avg = apply(l1,1,mean))
(nzo.avg = apply(nzo,1,mean))
(time.avg = apply(time,1,mean))

save(ans.avg, file=paste("ans", qv, ".rda", sep=""))
save(l1.avg, file=paste("l1", qv, ".rda", sep=""))
save(nzo.avg, file=paste("nzo", qv, ".rda", sep=""))
save(time.avg, file=paste("time", qv, ".rda", sep=""))

# q=0.5
#> (ans.avg = apply(ans,1,mean)*100)
#[1] 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333
#> (l1.avg = apply(l1,1,mean))
#[1] 0.4040006 0.4040006 0.4040006 0.3681103 0.3513791 0.2917211 0.2658054
#> (nzo.avg = apply(nzo,1,mean))
#[1]   3   4   4   6  29  85 135
#> (time.avg = apply(time,1,mean))
#[1] 24.229 23.433 31.616 60.689 66.951 43.982 33.337

# q=1
#> (ans.avg = apply(ans,1,mean)*100)
#[1] 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333
#> (l1.avg = apply(l1,1,mean))
#[1] 0.3681090 0.3681090 0.3681090 0.3681090 0.3354072 0.2784610 0.2537233
#> (nzo.avg = apply(nzo,1,mean))
#[1]   6   6   6   7  32 106 181
#> (time.avg = apply(time,1,mean))
#[1]  7.505  7.196  8.917 16.891 21.159 17.651 14.982


# q=2
#> (ans.avg = apply(ans,1,mean)*100)
#[1] 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333
#> (l1.avg = apply(l1,1,mean))
#[1] 0.3856286 0.3856286 0.3856286 0.3681011 0.3354001 0.2657988 0.2421860
#> (nzo.avg = apply(nzo,1,mean))
#[1]   7   7   7   9  35 128 209
#> (time.avg = apply(time,1,mean))
#[1]  4.636  4.442  4.948  7.367 11.457 10.200  9.513


# q=5
#> (ans.avg = apply(ans,1,mean)*100)
#[1] 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333 8.333333
#> (l1.avg = apply(l1,1,mean))
#[1] 0.3681142 0.3681142 0.3681142 0.3681142 0.3201670 0.2658082 0.2311864
#> (nzo.avg = apply(nzo,1,mean))
#[1]   7   7   7  11  43 134 247
#> (time.avg = apply(time,1,mean))
#[1]  4.506  4.383  4.495  5.272  8.085  9.259 10.422



