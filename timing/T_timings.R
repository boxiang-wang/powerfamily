rm(list=ls(all=TRUE))

if(.Platform$OS.type == "unix")
{
  file.loc = "/home/wang3660/Research/PF/timingtable/timing_0326/"
  setwd(file.loc)
  source("/home/wang3660/Research/PF/U_tool_server.R")
  require(gcdnet, Sys.getenv("R_LIBS_USER"))
}


if(.Platform$OS.type == "windows")
{
  file.loc = "D:/GitHub/powerfamily_output/all_timing/timing_0326/"
  setwd(file.loc)
  source("D:/GitHub/powerfamily/U_tool.R")
  require(gcdnet)
}

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

set.seed(1234)

if(!exists("nn")) nn=100
if(!exists("pp")) pp=5000
if(!exists("rho")) rho=0
if(!exists("total.indep")) total.indep=10
if(!exists("strong")) strong=F

seeds = sample(1:1000, size=1000)
FHT = FHTgen(n=nn, p=pp, rho=rho)
dat = FHT
x = dat$x
y = dat$y
seed.ind = 1

l2.list = c(0, 10^(-4), 10^(-2),1)
qv.list = c(0.5, 1, 2, 3, 5, 100)
avg.time.table = matrix(0, length(l2.list), (length(qv.list)+1))
if(strong == T) avg.time.tableF = matrix(0, length(l2.list), (length(qv.list)+1))


for(indp in 1:total.indep)
{
  print(paste(indp, nn, pp, rho, " th independent run.", sep=" "))
  time.table = matrix(NA, length(l2.list), (length(qv.list)+1))
  if(strong == T) time.tableF = matrix(NA, length(l2.list), (length(qv.list)+1))
  for(i in 1:length(l2.list)) # row
  {
    for(j in 1:(length(qv.list))) # column
    {
      
      l2 = l2.list[i]
      qv = qv.list[j]
      FHT = FHTgen(n=nn, p=pp, rho=rho)
      dat = FHT
      x = dat$x
      y = dat$y
      seed.ind = seed.ind + 1
      if(strong == T){
        start1 = Sys.time()
        m = GCDpower(x=x, y=y, lambda2=l2, qv=qv, method="power",
                     eps=1e-8, strong=F, standardize=T)
        stop1 = Sys.time()
        time.tableF[i,j] = difftime(stop1, start1, units="secs")
      }
      rm(m)
      start2 = Sys.time()
      m = GCDpower(x=x, y=y, lambda2=l2, qv=qv, method="power",
                   eps=1e-8, strong=T, standardize=T)
      stop2 = Sys.time()
      time.table[i,j] = difftime(stop2, start2, units="secs")
      
    }
    j = j + 1
	  start1 = Sys.time()
    m1 = gcdnet(x=x, y=y,lambda2=l2, method="logit",eps=1e-8, standardize=T)
    stop1 = Sys.time()
    time.table[i,j] = difftime(stop1, start1, units="secs")
    if(strong == T) time.tableF[i,j] = time.table[i,j]
  }
  # write.csv(time.table, file=paste("timetable_", indp, ".cvs", sep=""))
  #if (indp > 1)
  #{
  #  file.remove(paste("timetable_", indp - 1, ".rda", sep=""))
  #}
  avg.time.table = avg.time.table + time.table
  if(strong == T) avg.time.tableF = avg.time.tableF + time.tableF
}


avg.time.table = avg.time.table / total.indep
if(strong == T) avg.time.tableF = avg.time.tableF / total.indep
  
setwd(file.loc)

write.csv(avg.time.table, file=paste("avg.time", nn, pp, rho, ".csv", sep="_"))
save(avg.time.table, file=paste("avg.time", nn, pp, rho, ".rda", sep="_"))

if (strong == T){
  write.csv(avg.time.tableF, file=paste("avg.time", nn, pp, rho, "F.csv", sep="_"))
  save(avg.time.tableF, file=paste("avg.time", nn, pp, rho, "F.rda", sep="_"))
}

