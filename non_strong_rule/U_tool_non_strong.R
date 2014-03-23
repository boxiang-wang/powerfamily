setwd("D:\\GitHub\\powerfamily")
require(Matrix)
source("O_utilities.R")
source("M_p.GCDpower.R")
source("M_cv.GCDpower.R")
source("M_coef.GCDpower.R")
source("U_KKTcheckings.R")
source("M_FHTgen.R")


source("DrSVM/U_cv.DrSVM.fix2.R")
source("DrSVM/M_coef.DrSVM.R")
source("DrSVM/O_DrSVM_Fix2.R")
# dyn.load("M_powerfamilyNET.dll")

setwd("D:\\GitHub\\powerfamily\\non_strong_rule")



# dyn.unload("M_powerfamilyNET2.dll")
# shell("del M_powerfamilyNET2.dll M_powerfamilyNET2.o")
# shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET2.dll")
dyn.load("M_powerfamilyNET2.dll")

source("M_GCDpower_old.R")