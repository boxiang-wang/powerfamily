#rm(list=ls(all=TRUE))

simu.summary = function(nms, mm="mean"){
    
  if(!exists("nms")) nms='originalpath'
  print(nms) 
  
  if(.Platform$OS.type == "unix")
  {
    mainpath = "/home/wang3660/Research/PF/simu_accu/Mar29/"
    file.loc = paste(mainpath, nms, '/', sep="")
    setwd(file.loc)
    #source("/home/wang3660/Research/PF/U_tool_server.R")
    #require(gcdnet, Sys.getenv("R_LIBS_USER"))
  }
  
  
  if(.Platform$OS.type == "windows")
  {
    mainpath = "D:/GitHub/powerfamily_output/simu_accu/Mar29/"
    file.loc = paste(mainpath, nms, '/', sep="")
    if(nms %in% dir(mainpath) == FALSE) dir.create(file.loc) 
    setwd(file.loc)
    #source("D:/GitHub/powerfamily/U_tool.R")
    #require(gcdnet)
  }
  
  cv1st = c(6, 1:5)
  
  error_mm = matrix(NA, 6, 12)
  CI_mm = matrix(NA, 6, 12)
  
  se = function(x) sd(x)/sqrt(length(x))
  
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
  
  rownames(error_mm) = rownames(CI_mm) = c("q by CV", "q=0.5", "q=1", "q=2", "q=5", "q=100")
  
  error_mm = cbind(error_mm[, 1:4], NA, error_mm[, 5:8], NA, error_mm[, 9:12])
  colnames(error_mm)=c("l1", "l1_er", "enet", "enet_er",
                       NA, "l1", "l1_er", "enet", "enet_er", 
                       NA, "l1", "l1_er", "enet", "enet_er")
  
  CI_mm = cbind(CI_mm[, 1:4], NA, CI_mm[, 5:8], NA, CI_mm[, 9:12])
  colnames(CI_mm)=c("l1_C","l1_I","enet_C","enet_I",NA,"l1_C","l1_I","enet_C","enet_I",NA,
                    "l1_C","l1_I","enet_C","enet_I")
  info_mm=c(s_tr=info$s_tr, s_vl=info$s_vl, s_ts=info$s_ts, mu=info$mu,
            rho=info$rho)
  return(list(error_mm=error_mm, CI_mm=CI_mm, info_mm=info_mm))
}



options(scipen=10)
ss = simu.summary(55, "median")
ss$info_mm
ss$error_mm
ss$CI_mm




aa=function(x){
  for(j in c(2,4,7,9,12,14)){
     x[,j] = paste("(", x[,j], ")", sep="")
  }
  return(x)
}
ss$error_mm = round(ss$error_mm*100,2)
aa(ss$error_mm)


xtable(ss$error_mm, digits=2,sanitize.text.function = function(x){x})
print(xtable(aa(ss$error_mm), digits=2))


