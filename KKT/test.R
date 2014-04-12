rm(list=ls(all=T))

setwd("D:/GitHub/powerfamily/KKT")


dyn.unload("test2.dll")
shell("del test2.dll test2.o")
shell("Rcmd SHLIB test2.f90")

dyn.load("test2.dll")
a = 3.1
dd <- .Fortran("aa", a=a)


##################################################
rm(list=ls(all=T))


source("D:/GitHub/powerfamily/U_tool.R")
setwd("D:/GitHub/powerfamily/KKT")


qv = as.integer(2)
nn = as.integer(1000); pp = as.integer(5000)

set.seed(1234)
FHT = FHTgen(n=nn, p=pp, rho=0.8)
dat = FHT
x = dat$x; 
y = dat$y
lambda2 = as.double(1)

m = GCDpower(x=x,y=y,
             lambda=c(0.1,0.01),
             lambda2=lambda2, qv=2, method="power",eps=1e-4, standardize=F)
l = as.integer(length(m$lambda))

b0 = m$b0
betas = m$beta
rres = matrix(as.double(-1000), l, pp)
#k = 5
dyn.unload("test.dll")
shell("del test.dll test.o")
shell("Rcmd SHLIB test.f90")
dyn.load("test.dll")
thr = 1e-8
vio_count = as.integer(0)
dd <- .Fortran("KKT", q=qv, b0=b0, betas=as.matrix(betas),
               nn=nn, pp=pp, lam1=m$lambda, lam2=lambda2,
               l=l, x=x, y=y, thr=thr, vio_count=vio_count, quiet=T)
dd$vio_count/l/pp * 100


source("D:/GitHub/powerfamily/U_tool.R")
setwd("D:/GitHub/powerfamily/KKT")
KKTnp(b0, betas, y, x, c(0.1,0.01), lambda2, thr, loss = "power", qv=2)

r = y * (x %*% betas + matrix(rep(b0,nn),nn,length(m$lambda),byrow=T))
decib = qv/ (qv + 1)    
fdr =  decib ** (qv + 1.0)        

dlfun = function(r){
  ifelse(r > decib, (-1) * r ^ (-qv-1)*fdr, -1)
}

dl = apply(r, c(1,2), dlfun)

dly = dl * y
res = t(x) %*% (dl * y) / nn
res


lambda = m$lambda


count0 = 0
ctr0 = 0
ctrn0 = 0
for (l in 1:length(lambda)) {
  p = nrow(betas)
  for(j in 1:p)
  {
    if(betas[j,l]==0)
    {
      BB = abs(res[j,l]) - lambda[l]
      count0 = count0 + 1
      if (BB > thr) 
      {
        cat("violate at b = 0", BB, "\n")
        ctr0 <- ctr0 + 1
      }
    } else{
      AA = res[j,l] + lambda[l] * sign(betas[j,l]) + lambda2 * betas[j,l]
      if (abs(AA) >= thr)
      {
        cat("violate at b != 0", abs(AA), "\n")
        #print(betas[j,l])
        ctrn0 <- ctrn0 + 1
      }
      
    }
  }
}
ctr = ctrn0 + ctr0
ctr






## compare timing
load("D:\\GitHub\\powerfamily\\data\\prostate.rda")
x = prostate.x
y = prostate.y
dat = list(x = x, y = y)

eps = 1e-8
require(gcdnet)
system.time(m2 <- gcdnet(x=dat$x, y=dat$y, maxit=3e7,
                         lambda2=0, method="logit",eps=eps, standardize=T))[3]
m2$npass

require(glmnet)
system.time(m3 <- glmnet(x=dat$x, y=dat$y, maxit=3e7,
                         alpha=1, family="binomial",thresh=eps, standardize=T))[3]
m3$npass



##################################################
rm(list=ls(all=T))
source("D:/GitHub/powerfamily/U_tool.R")
setwd("D:/GitHub/powerfamily")

set.seed(1234)
FHT = FHTgen(n=100, p=3000, rho=0.5)
dat = FHT


load("D:\\GitHub\\powerfamily\\data\\arcene.rda")
x = tx
y = ty
dat = list(x = x, y = y)

dyn.unload("M_powerfamily.dll")
shell("del M_powerfamily.dll M_powerfamily.o")
shell("Rcmd SHLIB M_powerfamily.f90 O_auxiliary.f90 -o M_powerfamily.dll")
dyn.load("M_powerfamily.dll")

lambda2 = 1; qv = 2; eps = 1e-8; thr = 1e-3
#KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps, thr=thr)

system.time(m1 <- GCDpower(x=dat$x, y=dat$y,maxit=3e7,
                           lambda2=lambda2, qv=qv, method="power",eps=eps,
                           standardize=T, strong=T))[3]
m1$npass




source("D:/GitHub/powerfamily/non_strong_rule/U_tool_non_strong.R")
#KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps, thr=thr)

system.time(m2 <- GCDpower(x=dat$x, y=dat$y, maxit=3e7,
                  lambda2=lambda2, qv=qv, method="power",eps=eps, standardize=T))[3]
m2$npass


max(abs(m1$beta-m2$beta))









##################################################
rm(list=ls(all=T))

if(.Platform$OS.type == "unix")
{
  source("/home/wang3660/Research/PF/U_tool_server.R")
  setwd("/home/wang3660/Research/PF")
  load("/home/wang3660/Research/PF/data/arcene.rda")
}


if(.Platform$OS.type == "windows")
{
  source("D:/GitHub/powerfamily/U_tool.R")
  setwd("D:/GitHub/powerfamily")
  load("D:/GitHub/powerfamily/data/arcene.rda")
}


#load("D:\\GitHub\\powerfamily\\data\\arcene.rda")
#x = tx
#y = ty
#dat = list(x = x, y = y)

set.seed(1234)
FHT = FHTgen(n=5000, p=100, rho=0.5)
dat = FHT

lambda2 = 1; qv = 0.5; eps = 1e-8; thr = 1e-3
system.time(m1 <- GCDpower(x=dat$x, y=dat$y,maxit=3e7,
                           lambda2=lambda2, qv=qv, method="power",eps=eps,
                           standardize=T, strong=T))[3]
m1$npass


system.time(m2 <- GCDpower(x=dat$x, y=dat$y,maxit=3e7,
                           lambda2=lambda2, qv=qv, method="power",eps=eps,
                           standardize=T, strong=F))[3]
m2$npass


max(abs(m1$beta-m2$beta))






