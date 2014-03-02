###############################################################################
# Modified old function
rm(list=ls(all=T))

load("D:/GitHub/powerfamily/data/SPECTF.rda")
require(gcdnet)
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]


setwd("D:\\GitHub\\powerfamily\\DrSVM")


# Initialize the parameters in the function
nfolds = 5
delta = 2
pred.loss = "misclass"

N <- nrow(x)
y <- drop(y)

# Fit a model first, in order to get lambda sequence
GCDpower.object <- gcdnet(x, y, lambda = NULL, delta = delta)
lambda <- GCDpower.object$lambda

# record the number of non-zero coefficients for each lambda
nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)

# Start to initialize folder
set.seed(124)
foldid <- sample(rep(seq(nfolds), length = N)) 
outlist <- as.list(seq(nfolds))

# Within each folder, fit a model using training set
for (i in seq(nfolds)) {
  which <- foldid == i
  y_sub <- y[!which]
  outlist[[i]] <- gcdnet(x = x[!which, , drop = FALSE], 
                        y = y_sub, lambda = lambda, delta = delta, method="hhsvm")
}

# Compute the misclassification rate
predmat <- matrix(NA, length(y), length(lambda))
nlams <- double(nfolds)
for (i in seq(nfolds)) {
  which <- foldid == i
  fitobj <- outlist[[i]]
  preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
  nlami <- length(outlist[[i]]$lambda)
  predmat[which, seq(nlami)] <- preds
  nlams[i] <- nlami
}

# nfold rows and l column
# For each lambda, record misclassification rate for each folder
cvraw = (y != ifelse(predmat > 0, 1, -1))
outmat <- matrix(NA, nfolds, ncol(cvraw))
good <- matrix(0, nfolds, ncol(cvraw))
cvraw[is.infinite(cvraw)] <- NA
for (i in seq(nfolds)) {
  cvrawi <- cvraw[foldid == i, ]
  outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
  good[i, seq(nlams[i])] <- 1
}
cvraw = outmat

# For each folder, record the length of lambda
N <- apply(good, 2, sum)
cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
            cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
            nzero = nz, name = "Misclassification Error", 
            GCDpower.fit = GCDpower.object)
#lamin <- getmin(lambda, cvm, cvsd)

# Get the appropriate lambda
cvmin <- min(cvm)
idmin <- cvm <= cvmin
lambda.min <- max(lambda[idmin])
idmin <- match(lambda.min, lambda)
semin <- (cvm + cvsd)[idmin]
idmin <- cvm <= semin
lambda.1se <- max(lambda[idmin])
lamin <- list(lambda.min = lambda.min, lambda.1se = lambda.1se)
obj <- c(out, as.list(lamin))
class(obj) <- "cv.GCDpower"



#############################################################################
# Test new function cvs.GCDpower
rm(list=ls(all=T))

load("D:/GitHub/powerfamily/data/SPECTF.rda")
# require(gcdnet)
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]

setwd("D:\\GitHub\\powerfamily")
dyn.load("M_powerfamilyNET.dll")
source("M_GCDpower.R")
source("M_p.GCDpower.R")

setwd("D:\\GitHub\\powerfamily\\DrSVM")


cvs.GCDpower <- function(x, y, lambda = NULL, nfolds = 5, foldid, delta = 2, qv = 2, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  N <- nrow(x)
  y <- c(-1, 1)[as.factor(drop(y))]
  
  # Fit a model first, in order to get lambda sequence
  GCDpower.object <- GCDpower(x, y, lambda = lambda, delta = delta, qv = qv, 
                                ...)
  lambda <- GCDpower.object$lambda
  
  # record the number of non-zero coefficients for each lambda
  nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)
  
  # Start to initialize folder
  foldid <- sample(rep(seq(nfolds), length = N)) 
  outlist <- as.list(seq(nfolds))             
  
  # Setting some parameters
  predmat <- matrix(NA, length(y), length(lambda)) # record prediction
  nlams <- double(nfolds)                   # record lambda length for each folder
  
  # Treat each folder as test set, and observations out of this folder 
  #   as training set. Fit models using training sets, and record the prediction
  #   on test sets.
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- GCDpower(x = x[!which, , drop = FALSE], 
                           y = y[!which], lambda = lambda, delta = delta,
                                qv = qv, ...)
    preds <- predict(outlist[[i]], x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  
  # For each lambda, record misclassification rate for each folder
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  outmat <- matrix(NA, nfolds, ncol(cvraw))
  good <- matrix(0, nfolds, ncol(cvraw))
  cvraw[is.infinite(cvraw)] <- NA
  for (i in seq(nfolds)) {
    cvrawi <- cvraw[foldid == i, ]
    outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  cvraw = outmat
  
  N <- apply(good, 2, sum)  # For each folder, record the length of lambda
  
  # Compute the mean and the se of the mean for each lambda's cv
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
  # Find the largest lambda that can achieve the minimum cv
  #   and the largest lambda than can achieve the minimum cv plus 1 se of cv
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              nzero = nz, 
              lambda.min = lambda.min, lambda.1se = lambda.1se,
              name = "Misclassification Error", 
              GCDpower.fit = GCDpower.object)
  class(out) <- "cvs.GCDpower"
  out
}

set.seed(124)
m1 = cvs.GCDpower(x,y,method="hhsvm", delta=2)
m1$lambda.min

source("D:\\GitHub\\powerfamily\\M_cv.GCDpower.R")
set.seed(124)
m15 = cv.GCDpower(x,y,method="hhsvm", delta=2)
m15$lambda.min


set.seed(124)
require(gcdnet)
m2 = cv.gcdnet(x,y,method="hhsvm", delta=2, pred.loss="misclass")
m2$lambda.min



m1 = cv.gcdnet(x, y, delta = 2, method="hhsvm")
m2 = cv.GCDpower(x, y, delta = 2, method="hhsvm")


## for DrSVM
cvs.DrSVM_Fix2 <- function(x, y, lambda = NULL, lambda2 = 1, nfolds = 5, foldid, delta = 2, qv = 2, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  n <- nrow(x)
  y <- c(-1, 1)[as.factor(drop(y))]
  
  if(!is.null(lambda))
  {
    GCDlambda <- lambda
  }
  
  # Fit a model first, in order to get lambda sequence
  GCDpower.object <- GCDpower(x, y, lambda = lambda, delta = delta, qv = qv, 
                              ...)
  GCDlambda <- GCDpower.object$lambda
  
  # record the number of non-zero coefficients for each lambda
  nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)
  
  # Start to initialize folder
  foldid <- sample(rep(seq(nfolds), length = n)) 
  outlist <- as.list(seq(nfolds))             
  
  # Setting some parameters
  predmat <- matrix(NA, length(y), length(GCDpower)) # record prediction
  nlams <- double(nfolds)                   # record lambda length for each folder
  
  # Treat each folder as test set, and observations out of this folder 
  #   as training set. Fit models using training sets, and record the prediction
  #   on test sets.
  for (i in seq(nfolds)) {
    which <- foldid == i
    x_sub = x[!which, , drop = FALSE]
    outlist[[i]] <- DrSVM_Fix2(x = x_sub, y = y[!which], 
                              lambda2 = lambda2 * nrow(x_sub), ...)
    fitobj <- outlist[[i]]
    obj <- NULL
    obj$b0 <- (apply(fitobj$beta0, 1, mean))
    obj$beta <- t(fitobj$beta)
    obj$lambda <- fitobj$lambda1
    preds <- predict.GCDpower(obj, x[which, , drop = FALSE],
                              s=GCDlambda * nrow(x_sub), type="class")
    nlami <- length(outlist[[i]]$GCDlambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  
  # For each lambda, record misclassification rate for each folder
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  outmat <- matrix(NA, nfolds, ncol(cvraw))
  good <- matrix(0, nfolds, ncol(cvraw))
  cvraw[is.infinite(cvraw)] <- NA
  for (i in seq(nfolds)) {
    cvrawi <- cvraw[foldid == i, ]
    outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  cvraw = outmat
  
  # Fit DrSVM
  DrSVM.object <- DrSVM_Fix2(x, y, lambda2 = lambda2 * n, ...)
  
  # record the number of non-zero coefficients for each lambda
  nbeta <- t(DrSVM.object$beta)
  
  dimnames(nbeta) <- list(NULL, NULL)
  nzel <- function(x, which) if (any(x)) 
    which[x] else NULL
  t.beta <- abs(as.matrix(t(nbeta))) > 0
  
  nz <- (if (ncol(nbeta) == 1) 
    apply(t.beta, 2, nzel, seq(p)) else apply(t.beta, 1, nzel, seq(p)))
  nz <- sapply(nz, length)
  names(nz) <- paste("s", seq(nz), sep = "")
  
  
  # Create DrSVM fit
  newobj = NULL
  newobj$b0 <- (apply(DrSVM.object$beta0, 1, mean))
  newobj$beta <- nbeta
  newobj$lambda <- DrSVM.object$lambda1
  DrSVM.object = newobj
  
  N <- apply(good, 2, sum)  # For each folder, record the length of lambda
  
  # Compute the mean and the se of the mean for each lambda's cv
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
  # Find the largest lambda that can achieve the minimum cv
  #   and the largest lambda than can achieve the minimum cv plus 1 se of cv
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(GCDlambda[idmin])
  
  idmin <- match(lambda.min, GCDlambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(GCDlambda[idmin])
  
  out <- list(lambda = GCDlambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              nzero = nz, 
              lambda.min = lambda.min, lambda.1se = lambda.1se,
              name = "Misclassification Error", 
              DrSVM_Fix2.fit = DrSVM.object)
  class(out) <- "cvs.DrSVM_Fix2"
  out
}





###############################################################################
# Modify the function to apply for DrSVM_Fix2
rm(list=ls(all=T))

load("D:/GitHub/powerfamily/data/SPECTF.rda")
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]


setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Source files with tool functions.
source("O_utilities.R")

# Main program
source("M_GCDpower.R")

# Prediction, plot
source("M_p.GCDpower.R")
# KKT checking, CV
source("M_cv.GCDpower.R")
# coefficients
source("M_coef.GCDpower.R")
# KKT
source("U_KKTcheckings.R")
# Source file of data generator
source("M_FHTgen.R")
source("DrSVM/M_DrSVM_Fix2.R")


dyn.load("M_powerfamilyNET.dll")

# Initialize the parameters in the function
nfolds = 5
delta = 2
pred.loss = "misclass"
qv = 2 
lambda2= 1

lambda=NULL

# Initialize useful parameters
pred.loss = "misclass"


n <- nrow(x)
p <- ncol(x)
x <- as.matrix(x)
y <- c(-1, 1)[as.factor(drop(y))]

if(!is.null(lambda))
{
  GCDlambda <- lambda
}

# Fit a model first, in order to get lambda sequence
GCDpower.object <- GCDpower(x, y, lambda = lambda, delta = delta, qv = qv, standardize=F)
GCDlambda <- GCDpower.object$lambda
# record the number of non-zero coefficients for each lambda
nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)

# Start to initialize folder
foldid <- sample(rep(seq(nfolds), length = n)) 
outlist <- as.list(seq(nfolds))             

# Setting some parameters
predmat <- matrix(NA, length(y), length(GCDlambda)) # record prediction

nlams <- double(nfolds)                   # record lambda length for each folder

#   Treat each folder as test set, and observations out of this folder 
#   as training set. Fit models using training sets, and record the prediction
#   on test sets.
for (i in seq(nfolds)) {
  which <- foldid == i
  x_sub = x[!which, , drop = FALSE]
  outlist[[i]] <- DrSVM_Fix2(x = x_sub, y = y[!which], 
                             lambda2 = lambda2 * nrow(x_sub))
  fitobj <- outlist[[i]]
  obj <- NULL
  obj$b0 <- (apply(fitobj$beta0, 1, mean))
  obj$beta <- t(fitobj$beta)
  obj$lambda <- fitobj$lambda1
  preds <- predict.GCDpower(obj, x[which, , drop = FALSE],
                            s=GCDlambda * nrow(x_sub), type="class")
  nlami <- length(GCDlambda)
  predmat[which, seq(nlami)] <- preds
  nlams[i] <- nlami
}  
print(predmat)
# For each lambda, record misclassification rate for each folder
cvraw = (y != ifelse(predmat > 0, 1, -1))
outmat <- matrix(NA, nfolds, ncol(cvraw))
good <- matrix(0, nfolds, ncol(cvraw))
cvraw[is.infinite(cvraw)] <- NA
for (i in seq(nfolds)) {
  cvrawi <- cvraw[foldid == i, ]
  outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
  good[i, seq(nlams[i])] <- 1
}
cvraw = outmat

# Fit DrSVM
DrSVM.object <- DrSVM_Fix2(x, y, lambda2 = lambda2 * n, ...)

# record the number of non-zero coefficients for each lambda
nbeta <- t(DrSVM.object$beta)

dimnames(nbeta) <- list(NULL, NULL)
nzel <- function(x, which) if (any(x)) 
  which[x] else NULL
t.beta <- abs(as.matrix(t(nbeta))) > 0

nz <- (if (ncol(nbeta) == 1) 
  apply(t.beta, 2, nzel, seq(p)) else apply(t.beta, 1, nzel, seq(p)))
nz <- sapply(nz, length)
names(nz) <- paste("s", seq(nz), sep = "")


# Create DrSVM fit
newobj = NULL
newobj$b0 <- (apply(DrSVM.object$beta0, 1, mean))
newobj$beta <- nbeta
newobj$lambda <- DrSVM.object$lambda1
DrSVM.object = newobj

N <- apply(good, 2, sum)  # For each folder, record the length of lambda

# Compute the mean and the se of the mean for each lambda's cv
cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

# Find the largest lambda that can achieve the minimum cv
#   and the largest lambda than can achieve the minimum cv plus 1 se of cv
cvmin <- min(cvm)
idmin <- cvm <= cvmin
lambda.min <- max(GCDlambda[idmin])

idmin <- match(lambda.min, GCDlambda)
semin <- (cvm + cvsd)[idmin]
idmin <- cvm <= semin
lambda.1se <- max(GCDlambda[idmin])

out <- list(lambda = GCDlambda, cvm = cvm, cvsd = cvsd, 
            cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
            nzero = nz, 
            lambda.min = lambda.min, lambda.1se = lambda.1se,
            name = "Misclassification Error", 
            DrSVM_Fix2.fit = DrSVM.object)
class(out) <- "cvs.DrSVM_Fix2"
out

###############################################################################
# Test cvs.DrSVM_Fix2
rm(list=ls(all=T))
load("D:/GitHub/powerfamily/data/SPECTF.rda")
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]


setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Source files with tool functions.
source("O_utilities.R")

# Main program
source("M_GCDpower.R")

# Prediction, plot
source("M_p.GCDpower.R")
# KKT checking, CV
source("M_cv.GCDpower.R")
# coefficients
source("M_coef.GCDpower.R")
# KKT
source("U_KKTcheckings.R")
# Source file of data generator
source("M_FHTgen.R")
source("DrSVM/M_DrSVM_Fix2.R")


dyn.load("M_powerfamilyNET.dll")
load("D:/GitHub/powerfamily/data/SPECTF.rda")
require(gcdnet)
x = SPECTF.train[,-1]
x = as.matrix(x)
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]

#
load("data/colon.rda")
dat = NULL
dat$x = as.matrix(colon.x)
dat$y = c(-1,1)[as.factor(colon.y)]

x = as.matrix(dat$x); y = dat$y

qv = 5
n=nrow(x)

source("D:/GitHub/powerfamily/DrSVM/U_cv.DrSVM.fix2.R")
source("D:/GitHub/powerfamily/DrSVM/M_coef.DrSVM.R")
source("M_cv.GCDpower.R")

set.seed(124)
system.time(ccc1 <- cvs.DrSVM_Fix2(x, y, lambda2=0.1,scale=F))[3]
c1 = ccc1$coef.min

set.seed(124)
system.time(gg1 <- cv.GCDpower(x,y,lambda2=0.1, method="power", 
                               eps=1e-10, maxit=3e7, standardize=F))[3]
g1 = coef(gg1, s="lambda.min")

pre1 = sign(x%*%c1[-1]+matrix(rep(c1[1], n),nrow=n,byrow=T))
mean(y!=pre1)
pre2 = sign(x%*%g1[-1]+matrix(rep(g1[1], n),nrow=n,byrow=T))
mean(y!=pre2)

max(abs(c1[-1] - t(g1)[-1]))
pre1-pre2





dat = NULL
dat$x = matrix(rnorm(80*10),nrow=80)
dat$y = sign(dat$x[,1]+dat$x[,2]+dat$x[,3]+dat$x[,4])

system.time(ccc1 <- cvs.DrSVM_Fix2(dat$x, dat$y, lambda2=0.1,scale=F))[3]
#######################################################################
# test for classification rate table

system.time(tt <- DrSVM_Fix2(x, y, lambda2=1,scale=F))[3]
system.time(gg <- GCDpower(x,y,lambda2=1, method="power", eps=1e-10, 
                           maxit=3e7, standardize=F))[3]


lambda2 = 1e-4


lam = GCDpower(x, y, eps=1e-8, qv=2, delta=2,
               lambda2=1, method="power", standardize=F)$lambda
index = sample(1:nrow(x),as.integer(nrow(x)/5),replace=F) 
test_x = x[index,]
test_y = y[index]
train_x = x[-index,]
train_y = y[-index]

DrSVM_Fix2(train_x, train_y, lambda2 = lambda2 * nrow(train_x), smallmove=0.5)
