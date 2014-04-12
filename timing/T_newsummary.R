rm(list=ls(all=TRUE))

if(.Platform$OS.type == "unix")
{
  file.loc = "/home/wang3660/Research/PF/timingtable/timing_0323/"
  setwd(file.loc)
}


if(.Platform$OS.type == "windows")
{
  file.loc = "D:/GitHub/powerfamily_output/all_timing/timing_0326/"
  setwd(file.loc)
}

load("avg.time_100_5000_0_.rda")
m1 = avg.time.table
load("avg.time_100_5000_0.5_.rda")
m2 = avg.time.table
load("avg.time_100_5000_0.95_.rda")
m3 = avg.time.table

rbind(m1, m2, m3)

load("avg.time_100_5000_0_F.rda")
m1F = avg.time.tableF
load("avg.time_100_5000_0.5_F.rda")
m2F = avg.time.tableF
load("avg.time_100_5000_0.95_F.rda")
m3F = avg.time.tableF

rbind(m1F, m2F, m3F)

load("avg.time_5000_100_0_.rda")
m4 = avg.time.table
load("avg.time_5000_100_0.5_.rda")
m5 = avg.time.table
load("avg.time_5000_100_0.95_.rda")
m6 = avg.time.table

rbind(m4, m5, m6)

load("avg.time_5000_100_0_F.rda")
m4F = avg.time.tableF
load("avg.time_5000_100_0.5_F.rda")
m5F = avg.time.tableF
load("avg.time_5000_100_0.95_F.rda")
m6F = avg.time.tableF

rbind(m4F, m5F, m6F)



t1_no = (m1+m2+m3)/3
(a1 = apply(t1_no, 2, mean))

t1_F = (m1F+m2F+m3F)/3
(a1F = apply(t1_F, 2, mean))


t2_no = (m4+m5+m6)/3
(a2 = apply(t2_no, 2, mean))

t2_F = (m4F+m5F+m6F)/3
(a2F = apply(t2_F, 2, mean))


cc1 = mean(m1); cc2 = mean(m2); cc3 = mean(m3)
cc1F = mean(m1F); cc2F = mean(m2F); cc3F = mean(m3F)

cc4 = mean(m4); cc5 = mean(m5); cc6 = mean(m6)
cc4F = mean(m4F); cc5F = mean(m5F); cc6F = mean(m6F)

pdf("Comparison_T1_T2.pdf", 10, 10)
par(mfrow=c(1,2))

plot(1:6, a2F[-6], type="b", pch=1, ylim=c(0.96, 25), lty=2,
     xlab="", ylab="Timing (s)",xaxt='n', cex=1.2)
axis(side=1, at=1:6, labels=c("q=0.5", "DWD", "q=2", "q=3", "q=5", "logitnet"))

points(1:6, a1F[-6], type="b", pch=2, lty=2, cex=1.2)
points(1:6, a2[-6], type="b", pch=19, cex=1.2)
points(1:6, a1[-6], type="b", pch=17, cex=1.2)

leg.text = c("n=5000, p=100, non-strong", "n=5000, p=100, strong", "n=100, p=5000, non-strong", "n=100, p=5000, strong")
legend("topleft", pch=c(1,19,2,17), legend=leg.text, lty=c(2,1,2,1))

plot(1:3, c(cc4F, cc5F, cc6F), type="b", pch=1, ylim=c(2.5, 65), lty=2,
     xlab="", ylab="Timing (s)",xaxt='n')
axis(side=1, at=1:3, labels=c(expression(paste(rho, "= 0")), expression(paste(rho, "= 0.5")), expression(paste(rho, "= 0.95"))))
points(1:3, c(cc1F, cc2F, cc3F), type="b", pch=2, lty=2)

points(1:3, c(cc4, cc5, cc6), type="b", pch=19)
points(1:3, c(cc1, cc2, cc3), type="b", pch=17)
legend("topleft", pch=c(1,19,2,17), legend=leg.text, lty=c(2,1,2,1))
par(mfrow=c(1,1))
dev.off()



setwd("D:/Dropbox/Project/DWD/presentations/Presentation 04_03/image")
pdf("Comparison_T1_T2_noF.pdf", 10, 8)
par(mfrow=c(1,2))

plot(1:6, a2[-6], type="b", pch=19, ylim=c(0.96, 25), lty=2,
     xlab="", ylab="Timing (s)",xaxt='n', cex=1.5)
axis(side=1, at=1:6, labels=c("q=0.5", "DWD", "q=2", "q=3", "q=5", "logitnet"))
points(1:6, a1[-6], type="b", pch=17, cex=1.5, lty=3)


leg.text = c("n=5000, p=100", "n=100, p=5000")
legend("topleft", pch=c(19,17), legend=leg.text, lty=c(2,3))

plot(1:3, c(cc4, cc5, cc6), type="b", pch=19, ylim=c(2.5, 65), lty=2,
     xlab="", ylab="Timing (s)",xaxt='n', cex=1.5)
axis(side=1, at=1:3, labels=c(expression(paste(rho, "= 0")), expression(paste(rho, "= 0.5")), expression(paste(rho, "= 0.95"))))
points(1:3, c(cc1, cc2, cc3), type="b", pch=17, lty=3, cex=1.5)


legend("topleft", pch=c(19,17), legend=leg.text, lty=c(2,3))
par(mfrow=c(1,1))
dev.off()
