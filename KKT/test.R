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








##################################################
rm(list=ls(all=T))


source("D:/GitHub/powerfamily/U_tool.R")
setwd("D:/GitHub/powerfamily")



set.seed(1234)
FHT = FHTgen(n=80, p=95, rho=0.5)
dat = FHT



dyn.unload("M_powerfamilyNET.dll")
shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")
dyn.load("M_powerfamilyNET.dll")

lambda2 = 1; qv = 0.25; eps = 1e-6; thr = 1e-4
KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps, thr=thr)


dyn.unload("M_powerfamilyNET.dll")
shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
shell("Rcmd SHLIB M_powerfamilyNET2.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")
dyn.load("M_powerfamilyNET.dll")

lambda2 = 1; qv = 0.25; eps = 1e-6; thr = 1e-4
KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps, thr=thr)
