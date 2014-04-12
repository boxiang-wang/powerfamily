n = 50
blanc = 0.5
np = n * blanc
nm = n - np

p0 = 5
pp = 15

mu = 0.7
rho = 0.7

sina = 1

# The means of two classes.
mup = c(rep(mu, p0), rep(0, pp - p0))
mum = c(rep(-mu, p0), rep(0, pp - p0))

# The var-cov matrix of two classes.
if(sina == 1){
  sigma = diag(pp)
}
if(sina == 2){
  sigma = diag(pp)
  for(i in 1:p0) for(j in 1:p0) 
    sigma[i,j] = ifelse(i == j, 1, rho)
  
}
if(sina == 3){
  sigma = diag(pp)
  for(i in 1:p0) for(j in 1:p0) 
    sigma[i,j] = rho ^ abs(i-j)
}

eo = eigen(sigma, symmetric=TRUE)
sigma.sqrt = eo$vec %*% (diag(sqrt(eo$val)) %*% t(eo$vec))

plusgroup = 1:np
minusgroup = (np+1):n

x = matrix(NA, n, pp)
y = rep(NA, n)

Bayeserror.sum = 0
simu = 100000
for(tt in 1:simu){
  set.seed(tt)
  for(i in plusgroup){
    x[i,] = mup + sigma.sqrt %*% rnorm(pp)
    y[i] = 1
  }
  for( i in minusgroup){
    x[i,] = mum + sigma.sqrt %*% rnorm(pp)
    y[i] = -1 
  }
  
  
  Bayesrule = sign(apply(x[,1:p0], 1, sum))
  Bayeserror = mean(Bayesrule != y)
  Bayeserror.sum = Bayeserror.sum + Bayeserror
  
}

Bayeserror.sum/simu

