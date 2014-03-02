B = 30

dat = "arcene"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_27/", dat, sep=""))
q_arcene = rep(NA, B)
ans_arcene = matrix(NA, B, 8)
tim_arcene = matrix(NA, B, 8)

jj = c(1,2,4,5,6,8:18,20:23,25:30,103,107,119,180)
for(i in 1:B)
{
  
  load(paste("res_", dat, jj[i],".rda", sep=""))
  q_arcene[i] = prop.q
  ans_arcene[i, ] = ans.vec
  tim_arcene[i, ] = tim.vec
}

dat = "breast"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_27/", dat, sep=""))
q_breast = rep(NA, B)
ans_breast = matrix(NA, B, 8)
tim_breast = matrix(NA, B, 8)
for(i in 1:B)
{
  load(paste("res_", dat, i,".rda", sep=""))
  q_breast[i] = prop.q
  ans_breast[i, ] = ans.vec
  tim_breast[i, ] = tim.vec
}

dat = "colon"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_27/", dat, sep=""))
q_colon = rep(NA, B)
ans_colon = matrix(NA, B, 8)
tim_colon = matrix(NA, B, 8)
for(i in 1:B)
{
  load(paste("res_", dat, i,".rda", sep=""))
  q_colon[i] = prop.q
  ans_colon[i, ] = ans.vec
  tim_colon[i, ] = tim.vec
}

dat = "leuk"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_27/", dat, sep=""))
q_leuk = rep(NA, B)
ans_leuk = matrix(NA, B, 8)
tim_leuk = matrix(NA, B, 8)
for(i in 1:B)
{
  load(paste("res_", dat, i,".rda", sep=""))
  q_leuk[i] = prop.q
  ans_leuk[i, ] = ans.vec
  tim_leuk[i, ] = tim.vec
}

dat = "prostate"
setwd(paste("D:/GitHub/powerfamily/to server/crossvalidation/02_27/", dat, sep=""))
q_prostate = rep(NA, B)
ans_prostate = matrix(NA, B, 8)
tim_prostate = matrix(NA, B, 8)
for(i in 1:B)
{
  load(paste("res_", dat, i,".rda", sep=""))
  q_prostate[i] = prop.q
  ans_prostate[i, ] = ans.vec
  tim_prostate[i, ] = tim.vec
}


ww = c(1,3,4,5,7,8)

table(q_arcene)
p1 = ans_arcene[,ww]
p2 = tim_arcene[,ww]
A1 = apply(p1, 2, mean)
B1 = apply(p2, 2, mean)

table(q_breast)
p3 = ans_breast[,ww]
p4 = tim_breast[,ww]
A2 = apply(p3, 2, mean)
B2 = apply(p4, 2, mean)

table(q_colon)
p5 = ans_colon[,ww]
p6 = tim_colon[,ww]
A3 = apply(p5, 2, mean)
B3 = apply(p6, 2, mean)

table(q_leuk)
p7 = ans_leuk[,ww]
p8 = tim_leuk[,ww]
A4 = apply(p7, 2, mean)
B4 = apply(p8, 2, mean)

table(q_prostate)
p9 = ans_prostate[,ww]
p10 = tim_prostate[,ww]
A5 = apply(p9, 2, mean)
B5 = apply(p10, 2, mean)


A = rbind(A1,A2,A3,A4,A5)
B = rbind(B1,B2,B3,B4,B5)

tt = cbind(A[,1],B[,1],A[,2],B[,2],A[,3],B[,3],A[,4],B[,4],A[,5],B[,5],A[,6],B[,6])
require(xtable)
rownames(tt) = c("Arcene","Breast","Colon","Leukemia","Prostate")
xtable(tt, digits=3)


setwd("D:\\Dropbox\\Project\\DWD\\Draft\\images")
pdf("Plot_real_timings.pdf", 8, 8)

nf = layout(matrix(1:15, nrow=5,ncol=3, byrow=T))
#layout.show(nf)
o<-par(oma=c(3,0,3,0))

o<-par(mar=c(0,6,2,0))
barplot(table(q_arcene)/B, xaxt="n")
axis(3,at=1:5,labels=c(0.5,1,2,5,100))
mtext("Arcene", side=2, line=3)
mtext(expression(hat(q)), side=3, line=3)

o<-par(mar=c(0,2,2,0))
boxplot(p1, xaxt="n")
axis(3, at=1:6, labels=c(expression(hat(q)), "q=1", "q=2", "q=5", "SVM", "logis"))
mtext("Mis-Classification Rate", side=3, line=3)

o<-par(mar=c(0,2,2,4))
ylim=range(p2[,-5])*5/12
boxplot(p2[,-5], xaxt="n", yaxt="n", ylim=ylim)
axis(3, at=1:5, labels=c(expression(hat(q)), "q=1", "q=2", "q=5", "logis"))
axis(4, pretty(ylim,10))
mtext("Timings", side=3, line=3)



o<-par(mar=c(0,6,2,0))
barplot(table(q_breast)/B, xaxt="n")
mtext("Breast", side=2, line=3)
o<-par(mar=c(0,2,2,0))
boxplot(p3, xaxt="n")
o<-par(mar=c(0,2,2,4))
ylim=range(p4[,-5])/4
boxplot(p4[,-5], xaxt="n", yaxt="n", ylim=ylim)
axis(4, pretty(ylim,10))


o<-par(mar=c(0,6,2,0))
barplot(table(q_colon)/B, xaxt="n")
mtext("Colon", side=2, line=3)
o<-par(mar=c(0,2,2,0))
boxplot(p5, xaxt="n")
o<-par(mar=c(0,2,2,4))
ylim=range(p6[,-5])*2/3
boxplot(p6[,-5], xaxt="n", yaxt="n", ylim=ylim)
axis(4, pretty(ylim,10))


o<-par(mar=c(0,6,2,0))
barplot(table(q_leuk)/B, xaxt="n")
mtext("Leuk", side=2, line=3)
o<-par(mar=c(0,2,2,0))
boxplot(p7, xaxt="n")
o<-par(mar=c(0,2,2,4))
ylim=range(p8[,-5])/2
boxplot(p8[,-5], xaxt="n", yaxt="n", ylim=ylim)
axis(4, pretty(ylim,10))


o<-par(mar=c(0,6,2,0))
barplot(table(q_prostate)/B, xaxt="n")
mtext("Prostate", side=2, line=3)
o<-par(mar=c(0,2,2,0))
boxplot(p9, xaxt="n")
o<-par(mar=c(0,2,2,4))
ylim=range(p10[,-5])*2/3
boxplot(p10[,-5], xaxt="n", yaxt="n", ylim=ylim)
axis(4, pretty(ylim,10))

dev.off()


rrr = rbind(aaa1, aaa2, aaa3, aaa4, aaa5)
rrr = cbind(c(100,42,62,72,102), c(10000,22283,2000,7128,6033), rrr)
rownames(rrr) = c("Arcene","Breast","Colon","Leukemia","Prostate")
require(xtable)
xtable(rrr, digits=c(0,0,0,2,2,3,3,3))

rr = rbind(bb1, bb2, bb3, bb4, bb5)
rownames(rr) = c("Arcene","Breast","Colon","Leukemia","Prostate")
require(xtable)
xtable(rr, digits=3)

rrr1 = t(rr[,c(2,4,6,8,10)])
matplot(rrr1, type="l")