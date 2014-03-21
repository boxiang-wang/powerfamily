# Run GCDpower for the presentation
rm(list=ls(all=T))
setwd("/home/wang3660/Research/PF/timingtable/solnpathtime")
source("/home/wang3660/Research/PF/U_tool_server.R")


args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

dataname = "prostate"
if(dataname == "arcene")
{
  load("/home/wang3660/Research/PF/data/arcene.rda")
  x = tx
  y = ty
  y = c(-1,1)[as.factor(y)]
  l2 = 0.1
} 

if(dataname == "breast")
{
  load("/home/wang3660/Research/PF/data/breast.rda")
  y = c(-1,1)[as.factor(y)]
  l2 = 0.1
} 

if(dataname == "colon")
{
  load("/home/wang3660/Research/PF/data/colon.rda")
  x = colon.x
  y = colon.y
  y = c(-1,1)[as.factor(y)]
  l2 = 1
} 

if(dataname == "leuk")
{
  load("/home/wang3660/Research/PF/data/leuk.rda")
  y = c(-1,1)[as.factor(y)]
  l2 = 5
} 

if(dataname == "prostate")
{
  load("/home/wang3660/Research/PF/data/prostate.rda")
  x = prostate.x
  y = prostate.y
  y = c(-1,1)[as.factor(y)]
  l2 = 0.1
} 


m1 = GCDpower(x, y, method="power", qv=1, lambda2=1, standardize=F, eps=1e-08)
plot(m1)
title("q=1, lambda2=1", line = 3)

m2 = GCDpower(x, y, method="power", qv=2, lambda2=1, standardize=F, eps=1e-08)
plot(m2)
title("q=2, lambda2=1", line = 3)

m3 = GCDpower(x, y, method="power", qv=5, lambda2=1, standardize=F, eps=1e-08)
plot(m3)
title("q=5, lambda2=1", line = 3)


train_x =  x

train_m = apply(train_x, 2, mean)
std_train_x = t(apply(train_x, 1, function(x) x - train_m))  
train_sd = apply(std_train_x, 2, function(x) sqrt(x %*% x))
train_sd[train_sd==0] = 1
std_train_x = t(apply(std_train_x, 1, function(x) x / train_sd))  
train_x = std_train_x

x = train_x
x = prostate.x

tim = rep(NA, 20)
for(i in 1:20)
{
  tim[i] = system.time(fit1<-GCDpower(x, y, method="power", qv=1, 
                                    lambda2=1, standardize=F, 
                                    eps=1e-08))[3]
}
tim
median(tim)

tim = rep(NA, 20)
for(i in 1:20)
{
  tim[i] = system.time(fit1<-GCDpower(x, y, method="power", qv=1, 
                                      lambda2=1, standardize=T, 
                                      eps=1e-08))[3]
}
tim
median(tim)

x = train_x
fit1<-GCDpower(x, y, method="power", qv=1, 
               lambda2=1, standardize=T, 
               eps=1e-08)
fit2<-GCDpower(x, y, method="power", qv=1, 
               lambda2=1, standardize=F, 
               eps=1e-08)

dim(fit1$beta)
dim(fit2$beta)

start1=Sys.time()
md = DrSVM_Fix2(x, y, lambda2=5*nrow(x), eps=1e-07, scale=F)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

# cc = t(md$beta)
pdf("colonSVMpath.pdf", 10, 6)
cc = M_coef.powerfamily(newobj(md), s = m3$lambda*nrow(x))
s = apply(abs(cc[-1,]), 2, sum)
mm = cbind(t(cc[-1,]))
par(mar=c(5.1,5.1,4.1,2.1))
matplot(s, mm, type="l", cex.lab = 1.5, 
        lty=1, col = gray.colors(12, start = 0.05,
                                 end = 0.7, gamma = 2.2),
        
        xlim=c(0, max(s)), ylim=c(-0.15, 0.12),
        xlab=expression(paste("||",beta,"||",scriptscriptstyle(1))), 
        ylab=expression(beta), main="SVM lambda2=1 (Timing: 493s)")
dev.off()



pdf("colon_q_1.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(m1, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("q=1, lambda2=1 (Timing: 0.625s)", line = 3)
dev.off()

pdf("colon_q_2.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(m2, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("q=2, lambda2=1 (Timing: 0.485s)", line = 3)
dev.off()

pdf("colon_q_5.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(m3, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("q=5, lambda2=1 (Timing: 0.585s)", line = 3)
dev.off()



time1 = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  m1 = GCDpower(x, y, method="power", qv=1, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  time1[i]=difftime(stop1, start1, units="secs")
}
time2 = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  m2 = GCDpower(x, y, method="power", qv=2, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  time2[i]=difftime(stop1, start1, units="secs")
}
time5 = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  m3 = GCDpower(x, y, method="power", qv=5, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  time5[i]=difftime(stop1, start1, units="secs")
}
time1
mean(time1)

time2
mean(time2)

time5
mean(time5)


> time1
[1] 0.6202569 0.6235704 0.6104691 0.6179159 0.7017298 0.6128907 0.6128261
[8] 0.6138291 0.6210909 0.6176679
> mean(time1)
[1] 0.6252247
>
  > time2
[1] 0.4848092 0.4861219 0.4856305 0.4870214 0.4827800 0.4838407 0.4902256
[8] 0.4819033 0.4870102 0.4850717
> mean(time2)
[1] 0.4854414
>
  > time5
[1] 0.5763564 0.5754435 0.5780349 0.6535752 0.5728462 0.5770245 0.5820217
[8] 0.5758004 0.5789735 0.5844631
> mean(time5)
[1] 0.5854539



time05 = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  m05 = GCDpower(x, y, method="power", qv=0.5, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  time05[i]=difftime(stop1, start1, units="secs")
}
> time05
[1] 1.833414 1.835879 1.831553 1.896441 1.834627 1.833540 1.835989 1.831714
[9] 1.836842 1.832500
> mean(time05)
[1] 1.84025


time100 = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  m100 = GCDpower(x, y, method="power", qv=100, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  time100[i]=difftime(stop1, start1, units="secs")
}

> time100
[1] 2.857891 2.870060 2.847265 2.854826 2.860681 2.862548 2.856027 2.859933
[9] 2.859996 2.866434
> mean(time100)
[1] 2.859566


timeh = rep(NA, 10)
for(i in 1:10)
{
  
  start1=Sys.time()
  mh = GCDpower(x, y, method="hhsvm", delta=0.01, lambda2=1, standardize=F, eps=1e-07)
  stop1=Sys.time()
  timeh[i]=difftime(stop1, start1, units="secs")
}

> timeh
[1] 1.229869 1.266214 1.212815 1.224748 1.208614 1.206604 1.221035 1.216362
[9] 1.218320 1.219486
> mean(timeh)
[1] 1.222407




m05 = GCDpower(x, y, method="power", qv=0.5, lambda2=1, standardize=F, eps=1e-07)

m100 = GCDpower(x, y, method="power", qv=100, lambda2=1, standardize=F, eps=1e-07)

mh = GCDpower(x, y, method="hhsvm", delta=0.01, lambda2=1, standardize=F, eps=1e-07)


pdf("colon_q_05.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(m05, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("q=0.5, lambda2=1 (Timing: 1.840s)", line = 3)
dev.off()

pdf("colon_q_100.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(m100, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("q=100, lambda2=1 (Timing: 2.926s)", line = 3)
dev.off()

pdf("colon_hhsvm.pdf", 10, 6)
par(mar=c(5.1,5.1,4.1,2.1))
plot(mh, xlim=c(0, max(s)), ylim=c(-0.15, 0.12)
)
title("hhsvm, lambda2=1 (Timing: 1.222s)", line = 3)
dev.off()






cc1 = M_coef.powerfamily(newobj(md), s = (m100$lambda*nrow(x)))[-1,]
dimnames(cc1)=NULL
sum(abs(m100$beta - cc1))
sum(abs(mh$beta - cc1))
sum(abs(mb$beta - cc1))

mb = GCDpower(x, y, method="power", qv=120, lambda2=1, standardize=F, eps=1e-07)


sum(abs(mb$beta - mh$beta))
sum(abs(m3$beta - cc1))