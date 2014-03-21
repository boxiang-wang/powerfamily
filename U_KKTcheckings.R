# setwd("D:/GitHub/powerfamily/KKT")
# dyn.unload("KKT.dll")
# shell("del KKT.dll KKT.o")
# shell("Rcmd SHLIB KKT.f90")
dyn.load("D:/GitHub/powerfamily/KKT/KKT.dll")

# Usage:
# dd <- .Fortran("KKT", q=qv, b0=b0, betas=as.matrix(betas),
#               nn=nn, pp=pp, lam1=m$lambda, lam2=lambda2,
#               l=l, x=x, y=y, rres=rres, thr=thr, vio_count=vio_count, quiet=T)
# dd$vio_count/l/pp * 100

#################################################################################
############ Fit model and check KKT conditions ################
#################################################################################


KKTperctg = function(dat, lambda2, qv, eps, thr)
{
  if(eps > 1) eps = 10 ^ (-eps)
  if(thr > 1) thr = 10 ^ (-thr)
  
  m.temp = GCDpower(x=dat$x, y=dat$y,
                       lambda2=lambda2, qv=qv, method="power",eps=eps, standardize=F)
  
  vio_count = as.integer(0)
  dd = .Fortran("KKT", qv=qv, b0=m.temp$b0, betas=as.matrix(m.temp$beta),
                          nn=nrow(dat$x), pp=ncol(dat$x), lam1=m.temp$lambda, lam2=lambda2,
                          l=length(m.temp$lambda), x=dat$x, y=dat$y, 
                          thr=thr, vio_count=vio_count, quiet=T)
  dd$vio_count/dd$l/dd$pp*100 
}

#################################################################################
############ Summarize KKT condition checking tables ################
#################################################################################


KKTtb = function(dat, lambda2, qv, nm, eps.list=c(6:10), thr.list=c(2:5),
                 file.loc = "D:\\GitHub\\powerfamily_output\\KKT_check\\")
{
  perct.tb = matrix(NA, length(thr.list), length(eps.list))
  colnames(perct.tb) = paste("1e-", eps.list, sep="")
  rownames(perct.tb) = paste("1e-", thr.list, sep="")
  for(i in 1:length(eps.list))
  {
    for(j in 1:length(thr.list))
    {
      print(c(i, j))
      perct.tb[j, i] = KKTperctg(dat, lambda2=lambda2, qv=qv, 
                                 eps=eps.list[i], thr=thr.list[j])
    }
  }
  save("perct.tb", 
       file=paste(file.loc, nm, ".rda", sep="")) 
  return(perct.tb)
}
