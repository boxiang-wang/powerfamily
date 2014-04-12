#rm(list=ls(all=TRUE))

simu.summary = function(nms, mm="median"){
    
  if(!exists("nms")) nms='originalpath'
  print(nms) 
  
  if(.Platform$OS.type == "unix")
  {
    mainpath = "/home/wang3660/Research/PF/simu_accu_ada/Apr10/"
    file.loc = paste(mainpath, nms, '/', sep="")
    setwd(file.loc)
    #source("/home/wang3660/Research/PF/U_tool_server.R")
    #require(gcdnet, Sys.getenv("R_LIBS_USER"))
  }
  
  
  if(.Platform$OS.type == "windows")
  {
    mainpath = "D:/GitHub/powerfamily_output/simu_accu_ada/Apr10/"
    file.loc = paste(mainpath, nms, '/', sep="")
    if(nms %in% dir(mainpath) == FALSE) dir.create(file.loc) 
    setwd(file.loc)
    #source("D:/GitHub/powerfamily/U_tool.R")
    #require(gcdnet)
  }
  
  cv1st = c(6, 1:5)
  se = function(x) sd(x)/sqrt(length(x))
  
  # First, summarize the result without using adaptive methods.
  
  error_mm = matrix(NA, 6, 12)
  CI_mm = matrix(NA, 6, 12)
  
  load(file=paste(file.loc, "sina1/sina1allres.rda", sep=""))
  
  error_mm[,1] = apply(qverror_l1, 2, mm)[cv1st]
  error_mm[,2] = apply(qverror_l1, 2, se)[cv1st]
  CI_mm[,1] = apply(C_l1, 2, mm)[cv1st]
  CI_mm[,2] = apply(I_l1, 2, mm)[cv1st]
  
  error_mm[,3] = apply(qverror_elas, 2, mm)[cv1st]
  error_mm[,4] = apply(qverror_elas, 2, se)[cv1st]
  CI_mm[,3] = apply(C_elas, 2, mm)[cv1st]
  CI_mm[,4] = apply(I_elas, 2, mm)[cv1st]  
  
  load(file=paste(file.loc, "sina2/sina2allres.rda", sep=""))
  
  error_mm[,5] = apply(qverror_l1, 2, mm)[cv1st]
  error_mm[,6] = apply(qverror_l1, 2, se)[cv1st]
  CI_mm[,5] = apply(C_l1, 2, mm)[cv1st]
  CI_mm[,6] = apply(I_l1, 2, mm)[cv1st]
  
  error_mm[,7] = apply(qverror_elas, 2, mm)[cv1st]
  error_mm[,8] = apply(qverror_elas, 2, se)[cv1st]
  CI_mm[,7] = apply(C_elas, 2, mm)[cv1st]
  CI_mm[,8] = apply(I_elas, 2, mm)[cv1st]  
  
  load(file=paste(file.loc, "sina3/sina3allres.rda", sep=""))
  
  error_mm[,9] = apply(qverror_l1, 2, mm)[cv1st]
  error_mm[,10] = apply(qverror_l1, 2, se)[cv1st]
  CI_mm[,9] = apply(C_l1, 2, mm)[cv1st]
  CI_mm[,10] = apply(I_l1, 2, mm)[cv1st]
  
  error_mm[,11] = apply(qverror_elas, 2, mm)[cv1st]
  error_mm[,12] = apply(qverror_elas, 2, se)[cv1st]
  CI_mm[,11] = apply(C_elas, 2, mm)[cv1st]
  CI_mm[,12] = apply(I_elas, 2, mm)[cv1st]  
  
  # Then, summarize the result with using adaptive methods.
  
  error_mm_ada = matrix(NA, 6, 12)
  CI_mm_ada = matrix(NA, 6, 12)
  
  load(file=paste(file.loc, "sina1/sina1allres.rda", sep=""))
  
  error_mm_ada[,1] = apply(qverror_l1_ada, 2, mm)[cv1st]
  error_mm_ada[,2] = apply(qverror_l1_ada, 2, se)[cv1st]
  CI_mm_ada[,1] = apply(C_l1_ada, 2, mm)[cv1st]
  CI_mm_ada[,2] = apply(I_l1_ada, 2, mm)[cv1st]
  
  error_mm_ada[,3] = apply(qverror_elas_ada, 2, mm)[cv1st]
  error_mm_ada[,4] = apply(qverror_elas_ada, 2, se)[cv1st]
  CI_mm_ada[,3] = apply(C_elas_ada, 2, mm)[cv1st]
  CI_mm_ada[,4] = apply(I_elas_ada, 2, mm)[cv1st]  
  
  load(file=paste(file.loc, "sina2/sina2allres.rda", sep=""))
  
  error_mm_ada[,5] = apply(qverror_l1_ada, 2, mm)[cv1st]
  error_mm_ada[,6] = apply(qverror_l1_ada, 2, se)[cv1st]
  CI_mm_ada[,5] = apply(C_l1_ada, 2, mm)[cv1st]
  CI_mm_ada[,6] = apply(I_l1_ada, 2, mm)[cv1st]
  
  error_mm_ada[,7] = apply(qverror_elas_ada, 2, mm)[cv1st]
  error_mm_ada[,8] = apply(qverror_elas_ada, 2, se)[cv1st]
  CI_mm_ada[,7] = apply(C_elas_ada, 2, mm)[cv1st]
  CI_mm_ada[,8] = apply(I_elas_ada, 2, mm)[cv1st]  
  
  load(file=paste(file.loc, "sina3/sina3allres.rda", sep=""))
  
  error_mm_ada[,9] = apply(qverror_l1_ada, 2, mm)[cv1st]
  error_mm_ada[,10] = apply(qverror_l1_ada, 2, se)[cv1st]
  CI_mm_ada[,9] = apply(C_l1_ada, 2, mm)[cv1st]
  CI_mm_ada[,10] = apply(I_l1_ada, 2, mm)[cv1st]
  
  error_mm_ada[,11] = apply(qverror_elas_ada, 2, mm)[cv1st]
  error_mm_ada[,12] = apply(qverror_elas_ada, 2, se)[cv1st]
  CI_mm_ada[,11] = apply(C_elas_ada, 2, mm)[cv1st]
  CI_mm_ada[,12] = apply(I_elas_ada, 2, mm)[cv1st]  
  
  
  rownames(error_mm) = rownames(CI_mm) = c("$q$ by CV", "$q=0.5$", "$q=1$", "$q=2$", "$q=5$", "$q=100$")
  
  error_mm = cbind(error_mm[, 1:4], NA, error_mm[, 5:8], NA, error_mm[, 9:12])
  colnames(error_mm)=c("l1", "l1_er", "enet", "enet_er",
                       NA, "l1", "l1_er", "enet", "enet_er", 
                       NA, "l1", "l1_er", "enet", "enet_er")
  
  CI_mm = cbind(CI_mm[, 1:4], NA, CI_mm[, 5:8], NA, CI_mm[, 9:12])
  colnames(CI_mm)=c("l1_C","l1_I","enet_C","enet_I",NA,"l1_C","l1_I","enet_C","enet_I",NA,
                    "l1_C","l1_I","enet_C","enet_I")
  
  
  
  rownames(error_mm_ada) = rownames(CI_mm_ada) = c("$q$ by CV", "$q=0.5$", "$q=1$", "$q=2$", "$q=5$", "$q=100$")
  
  error_mm_ada = cbind(error_mm_ada[, 1:4], NA, error_mm_ada[, 5:8], NA, error_mm_ada[, 9:12])
  colnames(error_mm_ada)=c("l1", "l1_er", "enet", "enet_er",
                       NA, "l1", "l1_er", "enet", "enet_er", 
                       NA, "l1", "l1_er", "enet", "enet_er")
  
  CI_mm_ada = cbind(CI_mm_ada[, 1:4], NA, CI_mm_ada[, 5:8], NA, CI_mm_ada[, 9:12])
  colnames(CI_mm_ada)=c("l1_C","l1_I","enet_C","enet_I",NA,"l1_C","l1_I","enet_C","enet_I",NA,
                    "l1_C","l1_I","enet_C","enet_I")
  
  error = matrix(NA, 12, 11)
  rowind1 = c(1,3,5,7,9,11)
  rowind2 = c(2,4,6,8,10,12)
  error[rowind1,c(1,2,5,6,9,10)] = error_mm[,c(1,3,6,8,11,13)]
  error[rowind1,c(3,7,11)] = error_mm_ada[,c(3,8,13)]
  error[rowind2,c(1,2,5,6,9,10)] = error_mm[,c(2,4,7,9,12,14)]
  error[rowind2,c(3,7,11)] = error_mm_ada[,c(4,9,14)]
  
  rownames(error)  = c("$q$ by CV", NA, "$q=0.5$", NA, "$q=1$", NA, "$q=2$", NA, 
                          "$q=5$", NA, "$q=100$", NA)
  colnames(error)  = c("l1","enet","ada_enet",NA,"l1","enet","ada_enet",NA,
                       "l1","enet","ada_enet")
  
  CI = matrix(NA, 6, 20)
  CI = cbind(CI_mm[, 1:4],  CI_mm_ada[, 3:4],NA, 
                    CI_mm[, 6:9],  CI_mm_ada[, 8:9],NA, 
                    CI_mm[, 11:14],  CI_mm_ada[, 13:14])
  colnames(CI)=c("l1_C","l1_I","enet_C","enet_I","ada_enet_C","ada_enet_I",NA,
                        "l1_C","l1_I","enet_C","enet_I","ada_enet_C","ada_enet_I",NA,
                        "l1_C","l1_I","enet_C","enet_I","ada_enet_C","ada_enet_I")
  
  info_mm=c(s_tr=info$s_tr, s_vl=info$s_vl, s_ts=info$s_ts, mu=info$mu,
            rho=info$rho)
  return(list(error_mm=error_mm, CI_mm=CI_mm, 
              error_mm_ada=error_mm_ada, CI_mm_ada=CI_mm_ada, 
              info_mm=info_mm,
              error=error,CI=CI
              ))
}


options(scipen=10)
ss = simu.summary(17, "median")
ss$info_mm
round(ss$error*100,4)[c(1,3,5,7,9,11),]
ss$CI


# 17, 22, 23 are good.
# 25 is very good
# 27 is very good

ss$CI_mm
ss$error_mm_ada
ss$CI_mm_ada
ss$error_mm

aa=function(x){
  x = x*100
  cin = c(1:3,5:7,9:11)
  #x=as.data.frame(x)
  for(j in c(1:6)*2){
     tp = as.numeric(x[j,cin])
     x[j,cin] = paste0("(", sprintf("%.2f",tp), ")", sep="")
  }
  for(j in c(1,3,5,7,9,11)){
     tp = as.numeric(x[j,cin])
     x[j,cin] = sprintf("%.2f",tp)
  }
  return(x)
}

kk = aa(ss$error)
kk2 = cbind(rownames(kk), kk)
colnames(kk2) = c(NA, colnames(kk))
rownames(kk2) = NULL

require(xtable)
print(xtable(kk2),
      digits=2, booktabs=T, sanitize.text.function = function(x){x}, include.rownames=F)

ssCI2 = ss$CI
for(ii in 1:6){
  ssCI2[ii,] = paste(round(ss$CI[ii,],1))
}


print(xtable(ssCI2),
      booktabs=T, sanitize.text.function = function(x){x}, include.rownames=F)
