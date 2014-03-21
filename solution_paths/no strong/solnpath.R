rm(list=ls(all=T))

# if l2_slct = 1, we use lambda2 = 1 all the time.
# if l2_slct = 2, we use lambda2 from the cross validation.

# if x_slct = 1, we use the full data. 
# if x_slct = 2, we use training data.

lambda2seq = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10)
qvseq = c(0.5, 1, 2, 5, 100)
totaltime = 10

dat = "colon"
rs = 4
oneplot = TRUE
x_slct = 1
l2_slct = 1
fix_lam = 0
eps = 1e-8

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if(.Platform$OS.type == "unix")
{
  setwd(paste("/home/wang3660/Research/PF/timingtable/solnpathtime/no_strong/", dat, sep=""))
  source("/home/wang3660/Research/PF/U_tool_server.R")
  datapath = "/home/wang3660/Research/PF/data/"
  require(gcdnet, Sys.getenv("R_LIBS_USER"))
}


if(.Platform$OS.type == "windows")
{
  setwd(paste("D:/GitHub/powerfamily_output/solutionpath/no_strong/", dat, sep=""))
  source("D:/GitHub/powerfamily/U_tool.R")
  datapath = "D:/GitHub/powerfamily/data/"
  require(gcdnet)
}

load(paste(datapath, dat, ".rda", sep=""))

if(dat == "arcene")
{
  x = tx
  y = ty
  y = c(-1,1)[as.factor(y)]
} 

if(dat == "breast")
{
  y = c(-1,1)[as.factor(y)]
} 

if(dat == "colon")
{
  x = colon.x
  y = colon.y
  y = c(-1,1)[as.factor(y)]
} 

if(dat == "leuk")
{
  y = c(-1,1)[as.factor(y)]
} 

if(dat == "prostate")
{
  x = prostate.x
  y = prostate.y
  y = c(-1,1)[as.factor(y)]
} 

set.seed(rs)
x = as.matrix(x)
y = c(-1,1)[as.factor(y)]

set.seed(rs)
index = sample(1:nrow(x),as.integer(nrow(x)/3),replace=F) 
test_x = x[index,]
test_y = y[index]
train_x = x[-index,]
train_y = y[-index]


if.std = 1
if(if.std == 1)
{
  train_m = apply(train_x, 2, mean)
  std_train_x = t(apply(train_x, 1, function(x) x - train_m))  
  train_sd = apply(std_train_x, 2, function(x) sqrt(x %*% x / length(x)))
  train_sd[train_sd==0] = 1
  std_train_x = t(apply(std_train_x, 1, function(x) x / train_sd))  
  
  std_test_x = t(apply(test_x, 1, function(x) x - train_m))  
  std_test_x = t(apply(std_test_x, 1, function(x) x / train_sd)) 
  
  train_x = std_train_x
  test_x = std_test_x
}
if(if.std == 1)
{
  x_m = apply(x, 2, mean)
  std_x = t(apply(x, 1, function(x) x - x_m))  
  x_sd = apply(std_x, 2, function(x) sqrt(x %*% x / length(x)))
  x_sd[x_sd==0] = 1
  std_x = t(apply(std_x, 1, function(x) x / x_sd))  
  
}

###################################TUNED##########################################
# if x_slct = 1, we use the full data. 
# if x_slct = 2, we use training data.

if(x_slct == 1){
  working_x = std_x
  working_y = y
}
if(x_slct == 2){
  working_x = train_x
  working_y = train_y
}


##############
# if l2_slct = 1, we use lambda2 = fix.lam all the time.
# if l2_slct = 2, we use lambda2 from the cross validation.

if(l2_slct == 2)
{
  # First, we need to find the best lambda_2 for each model
  nrep = length(qvseq)
  ans = matrix(0, length(lambda2seq), nrep)
  j = 1
  for(qv in qvseq)
  {
    i = 1
    print(i)
    for(lambda2 in lambda2seq){
      set.seed(rs)
      cv<-cv.GCDpower(working_x, working_y, eps=eps, qv=qv,
                      lambda2=lambda2, method="power",
                      standardize = FALSE,
                      pred.loss="misclass", nfolds=5)
      cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
      ans[i,j] = cvm
      i = i + 1
    }
    j = j + 1
  }
  
  # These values are the best values for each fixed q
  bestl2_power = apply(ans, 2, function(x) length(lambda2seq) - which.min(rev(x)) + 1 )
  cvl2 = lambda2seq[bestl2_power]
  
  
  # Tune lambda2 for logistic regression
  ans1 = rep(NA, (length(lambda2seq)))
  start2 = Sys.time()
  i = 1
  for(lambda2 in lambda2seq){
    set.seed(rs)
    fittime = system.time(cv<-cv.gcdnet(working_x, working_y, eps=eps, 
                                        lambda2=lambda2, method="logit",
                                        standardize = FALSE,
                                        pred.loss="misclass", nfolds=5))[3]
    cvm = cv$cvm[which(cv$lambda==cv$lambda.1se)]
    ans1[i] = cvm
    i = i + 1
  }
  
  bestl2_logis = length(lambda2seq) - which.min(rev(ans1)) + 1
  cvl2_logis = lambda2seq[bestl2_logis]
  
}

if(l2_slct == 1)
{
  cvl2_working = rep(fix_lam, length(qvseq))
  cvl2_working_logis = fix_lam
}

if(l2_slct == 2)
{
  cvl2_working = cvl2
  cvl2_working_logis = cvl2_logis
}

#############################################
# Then we get the time
lf = c(qvseq, "SVM", "logis_gcd")
N = length(lf)

tim = matrix(NA, N, totaltime)
for(l in 1:length(qvseq)){
  for(i in 1:totaltime){
    tim[l, i] = system.time(fit<-GCDpower(working_x, working_y, method="power", qv=qvseq[l], 
                                lambda2=cvl2_working[l], standardize=F, eps=eps))[3]
  }
}

# For SVM
#for(i in 1:totaltime){
#  tim[6, i] 
if(fix_lam != 0)
{
  tim[6,] = system.time(fit_D<-DrSVM_Fix2(working_x, working_y, eps=eps, 
                                          lambda2=cvl2_working[length(qvseq)]*nrow(working_x),
                                          scale=F))[3]
}

#  print(i)
#}

# For logistic (gcdnet) 
for(i in 1:totaltime){
  tim[7, i] = system.time(fit_logis_gcdnet<-gcdnet(working_x, working_y, eps=1e-8, 
                                        lambda2=cvl2_working_logis, method="logit", 
                                        standardize=F))[3]
}
avg.time = apply(tim, 1, mean)



if(fix_lam == 0)
{
  require(glmnet)
  tim_glmnet = rep(NA, totaltime)
  for(i in 1:totaltime){
    tim_glmnet[i] = system.time(fit_logis_glmnet<-
    glmnet(working_x, working_y, thresh=eps, family="binomial",
           alpha=1, standardize=F,lambda=fit_logis_gcdnet$lambda))[3]
  }
}
####################################PLOT##########################################
for(oneplot in c(T,F))
{
  # This is the average time
  tim_plot = avg.time
  
  # Decide what is plotted
  plot_x = working_x
  plot_y = working_y
  
  # Size
  width = 7
  height = 4.2
  
  if(l2_slct == 1){ #we use lambda2 = 1 all the time.
    if(oneplot == T){
      pdf(paste(dat, "_lam_", fix_lam, ".pdf", sep=""), width, height)
    }
  }
  if(l2_slct == 2){
    if(oneplot == T){
      pdf(paste(dat, ".pdf", sep=""), width, height)
    }
  }
  
  
  # Plot models
  # First decide the plotting range by SVM
  if(fix_lam != 0){
    if(!exists('fit_D')){
      fit_D = DrSVM_Fix2(plot_x, plot_y, eps=eps, 
                         lambda2=cvl2_working[length(qvseq)]*nrow(plot_x),
                         scale=F)
    }
    cc = M_coef.powerfamily(newobj(fit_D), s = fit$lambda*nrow(plot_x))
    s = apply(abs(cc[-1,]), 2, sum)
    xlim = c(0, max(s)) * 0.9
    mm = cbind(t(cc[-1,]))
    r1 = max(mm) * 1.1
    r2 = min(mm) * 1.1
    ylim = c(r2, r1)
  }
  if(fix_lam == 0){
    m100 = mfit = GCDpower(plot_x, plot_y, method="power", qv=100, lambda2=fix_lam, 
                           standardize=F, eps=eps)
    cc = M_coef.powerfamily(mfit, s = mfit$lambda)
    s = apply(abs(cc[-1,]), 2, sum)
    xlim = c(0, max(s)) * 0.9
    mm = t(cc[-1,])
    r1 = max(mm) * 1.1
    r2 = min(mm) * 1.1
    ylim = c(r2, r1) 
  }
  
  
  i = 1
  m05 = mfit = GCDpower(plot_x, plot_y, method="power", qv=qvseq[i], lambda2=cvl2_working[i], 
                        standardize=F, eps=eps)
  if(oneplot == F) pdf(paste(dat, "_q_", qvseq[i],
                             "_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit
       , xlim = xlim, ylim = ylim
  )
  title(paste("q = ", qvseq[i], ", lambda2 = ", 
              cvl2_working[i], " (Timing: ", round(tim_plot[i], 3), " s)", sep=""), line = 3)
  if(oneplot == F) dev.off()
  
  i = 2
  m1 = mfit = GCDpower(plot_x, plot_y, method="power", qv=qvseq[i], lambda2=cvl2_working[i], 
                       standardize=F, eps=eps)
  if(oneplot == F) pdf(paste(dat, "_q_", qvseq[i],
                             "_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit
       , xlim = xlim, ylim = ylim
  )
  title(paste("q = ", qvseq[i], ", lambda2 = ", 
              cvl2_working[i], " (Timing: ", round(tim_plot[i], 3), " s)", sep=""), line = 3)
  if(oneplot == F) dev.off()
  
  i = 3
  m2 = mfit = GCDpower(plot_x, plot_y, method="power", qv=qvseq[i], lambda2=cvl2_working[i], 
                       standardize=F, eps=eps)
  if(oneplot == F) pdf(paste(dat, "_q_", qvseq[i],
                             "_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit
       , xlim = xlim, ylim = ylim
  )
  title(paste("q = ", qvseq[i], ", lambda2 = ", 
              cvl2_working[i], " (Timing: ", round(tim_plot[i], 3), " s)", sep=""), line = 3)
  if(oneplot == F) dev.off()
  
  
  i = 4
  m5 = mfit = GCDpower(plot_x, plot_y, method="power", qv=qvseq[i], lambda2=cvl2_working[i], 
                       standardize=F, eps=eps)
  if(oneplot == F) pdf(paste(dat, "_q_", qvseq[i],
                             "_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit
       , xlim = xlim, ylim = ylim
  )
  title(paste("q = ", qvseq[i], ", lambda2 = ", 
              cvl2_working[i], " (Timing: ", round(tim_plot[i], 3), " s)", sep=""), line = 3)
  if(oneplot == F) dev.off()
  
    
  i = 5
  m100 = mfit = GCDpower(plot_x, plot_y, method="power", qv=qvseq[i], lambda2=cvl2_working[i], 
                         standardize=F, eps=eps)
  if(oneplot == F) pdf(paste(dat, "_q_", qvseq[i],
                             "_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit
       , xlim = xlim, ylim = ylim
  )
  title(paste("q = ", qvseq[i], ", lambda2 = ", 
              cvl2_working[i], " (Timing: ", round(tim_plot[i], 3), " s)", sep=""), line = 3)
  if(oneplot == F) dev.off()
  
  
  if(fix_lam != 0){
    # SVM
    mfit = fit_D
    if(oneplot == F) pdf(paste(dat, "_SVM_lam_", cvl2_working[i], ".pdf", sep=""), width, height)
    par(mar=c(5.1,5.1,4.1,2.1))
    matplot(s, mm, type="l", cex.lab = 1, pch = 500,
            lty=1, col = gray.colors(12, start = 0.05, end = 0.7, gamma = 2.2),
            xlim = xlim, ylim = ylim, 
            xlab="L1 Norm", 
            ylab="Coefficients")
    title(paste("SVM, lambda2 = ", 
                cvl2_working[length(qvseq)], " (Timing: ", round(tim_plot[length(qvseq)+1], 3), " s)",
                sep=""), line = 3)
    if(oneplot == F) dev.off()
  }
  
  # logistic
  mfit = fit_logis_gcdnet = gcdnet(plot_x, plot_y, eps=eps, 
                                   lambda2=cvl2_working_logis, method="logit", 
                                   standardize=F)
  if(oneplot == F) pdf(paste(dat, "_logistic(gcdnet)_lam_", 
                             cvl2_working_logis, ".pdf", sep=""), width, height)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(mfit, xlim = xlim, ylim = ylim)
  title(paste("Logistic, lambda2 = ", 
              cvl2_working_logis, " (Timing: ", round(tim_plot[N], 3), " s)",
              sep=""), line = 3)
  if(oneplot == F) dev.off()
  
  if(fix_lam == 0)
  {
    # logistic (glmnet)
    mfit = fit_logis_glmnet = glmnet(plot_x, plot_y, thresh=eps, family="binomial",
                                     alpha=1, standardize=F,
                                     lambda=fit_logis_gcdnet$lambda)
    if(oneplot == F) pdf(paste(dat, "_logistic(glmnet)_lam_", 
                               cvl2_working_logis, ".pdf", sep=""), width, height)
    par(mar=c(5.1,5.1,4.1,2.1))
    plot(mfit, xlim = xlim, ylim = ylim)
    title(paste("Logistic (glmnet), lambda2 = ", 
                cvl2_working_logis, " (Timing: ", round(mean(tim_glmnet), 3), " s)",
                sep=""), line = 3)
    if(oneplot == F) dev.off()
  }
  
  
  if(oneplot == T)
  {
    dev.off()
  }
}




if(1 == 2){
  require(glmnet)
  mfit = fit_logis_glmnet = glmnet(plot_x, plot_y, thresh=eps, family="binomial",
                                   alpha=1, standardize=F,lambda=fit_logis_gcdnet$lambda)
  plot(mfit, ylim = c(-2,2))
  
  mfit = fit_logis_gcdnet = gcdnet(plot_x, plot_y, eps=eps, 
                                   lambda2=0, method="logit", 
                                   standardize=F)
  plot(mfit, ylim = c(-2,2))
}