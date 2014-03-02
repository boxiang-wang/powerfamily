source("D:/GitHub/powerfamily/DrSVM/O_DrSVM_Fix2.r")

load("D:/GitHub/powerfamily/data/SPECTF.rda")
x = as.matrix(SPECTF.train[,-1])
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]

index = sample(1:nrow(x),as.integer(nrow(x)/5),replace=F) 
test_x = x[index,]
test_y = y[index]
train_x = x[-index,]
train_y = y[-index]


m = DrSVM_Fix2(train_x, train_y, lambda2=1)

load("/home/wang3660/Research/PF/data/arcene.rda")
source("/home/wang3660/Research/PF/classrate/DrSVM/O_DrSVM_Fix2.R")
dat = NULL
dat$x = tx
dat$y = ty
dat$y = c(-1, 1)[as.factor(dat$y)]

m = DrSVM_Fix2(dat$x, dat$y, lambda2=0.1)