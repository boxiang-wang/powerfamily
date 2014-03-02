cvs.GCDpower <- function(x, y, lambda = NULL, nfolds = 5, foldid, delta = 2, qv = 2, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  n <- nrow(x)
  y <- c(-1, 1)[as.factor(drop(y))]
  
  # Fit a model first, in order to get lambda sequence
  GCDpower.object <- GCDpower(x, y, lambda = lambda, delta = delta, qv = qv, 
                                ...)
  lambda <- GCDpower.object$lambda
  
  # record the number of non-zero coefficients for each lambda
  nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)
  
  # Start to initialize folder
  foldid <- sample(rep(seq(nfolds), length = n)) 
  outlist <- as.list(seq(nfolds))             
  
  # Setting some parameters
  predmat <- matrix(NA, length(y), length(lambda)) # record prediction
  nlams <- double(nfolds)                   # record lambda length for each folder
  
  # Treat each folder as test set, and observations out of this folder 
  #   as training set. Fit models using training sets, and record the prediction
  #   on test sets.
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- GCDpower(x = x[!which, , drop = FALSE], 
                           y = y[!which], lambda = lambda, delta = delta,
                                qv = qv, ...)
    preds <- predict(outlist[[i]], x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  
  # For each lambda, record misclassification rate for each folder
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  outmat <- matrix(NA, nfolds, ncol(cvraw))
  good <- matrix(0, nfolds, ncol(cvraw))
  cvraw[is.infinite(cvraw)] <- NA
  for (i in seq(nfolds)) {
    cvrawi <- cvraw[foldid == i, ]
    outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  cvraw = outmat
  
  N <- apply(good, 2, sum)  # For each folder, record the length of lambda
  
  # Compute the mean and the se of the mean for each lambda's cv
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
  # Find the largest lambda that can achieve the minimum cv
  #   and the largest lambda than can achieve the minimum cv plus 1 se of cv
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              nzero = nz, 
              lambda.min = lambda.min, lambda.1se = lambda.1se,
              name = "Misclassification Error", 
              GCDpower.fit = GCDpower.object)
  class(out) <- "cvs.GCDpower"
  out
}

cvs.DrSVM_Fix2 <- function(x, y, lambda2 = 1, lambda = NULL, nfolds = 5, foldid, delta = 2, qv = 2, scale=F, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  n <- nrow(x)
  p <- ncol(x)
  x <- as.matrix(x)
  y <- c(-1, 1)[as.factor(drop(y))]
  
  if(!is.null(lambda))
  {
    GCDlambda <- lambda
  }
  
  # Fit a model first, in order to get lambda sequence
  GCDpower.object <- GCDpower(x, y, lambda = lambda, delta = delta, qv = qv, standardize=scale)
  GCDlambda <- GCDpower.object$lambda
  # Fit DrSVM
  DrSVM.object <- DrSVM_Fix2(x, y, lambda2 = lambda2 * n, ...)
  signal = TRUE
  while(signal == TRUE)
  {
    signal = FALSE
	# Start to initialize folder
    foldid <- sample(rep(seq(nfolds), length = n)) 
    outlist <- as.list(seq(nfolds))             
  
    # Setting some parameters
    predmat <- matrix(NA, length(y), length(GCDlambda)) # record prediction

    nlams <- double(nfolds)                   # record lambda length for each folder  

    #   Treat each folder as test set, and observations out of this folder 
    #   as training set. Fit models using training sets, and record the prediction
    #   on test sets.
    for (i in seq(nfolds)) {
      which <- foldid == i
      x_sub = x[!which, , drop = FALSE]
      tempmodel <- DrSVM_Fix2(x = x_sub, y = y[!which], 
                                 lambda2 = lambda2 * nrow(x_sub), ...)
      if(class(tempmodel) != "list")
  	  {
	    signal = TRUE
		print('resample3')
		break
	  }
	  if(exists('tempmodel$signal'))
	  {
	    if(tempmodel$signal == TRUE)
	    {
	      signal = TRUE
		  print('resample3')
		  break
	    }
	  }
	  outlist[[i]] <- tempmodel
	  fitobj <- outlist[[i]]
      obj <- NULL
      obj$b0 <- (apply(fitobj$beta0, 1, mean))
      obj$beta <- t(fitobj$beta)
      obj$lambda <- fitobj$lambda1
      preds <- predict.GCDpower(obj, x[which, , drop = FALSE],
                                s=GCDlambda * nrow(x_sub), type="class")
      predmat[which, seq(GCDlambda)] <- preds
      nlams[i] <- length(GCDlambda)
    }
  } 
  # For each lambda, record misclassification rate for each folder
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  outmat <- matrix(NA, nfolds, ncol(cvraw))
  good <- matrix(0, nfolds, ncol(cvraw))
  cvraw[is.infinite(cvraw)] <- NA
  for (i in seq(nfolds)) {
    cvrawi <- cvraw[foldid == i, ]
    outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  cvraw = outmat
  
  
  # Create DrSVM fit
  newobj = NULL
  newobj$b0 <- (apply(DrSVM.object$beta0, 1, mean))
  newobj$beta <- t(DrSVM.object$beta)
  newobj$lambda <- DrSVM.object$lambda1
  
  N <- apply(good, 2, sum)  # For each folder, record the length of lambda
  
  # Compute the mean and the se of the mean for each lambda's cv
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
  # Find the largest lambda that can achieve the minimum cv
  #   and the largest lambda than can achieve the minimum cv plus 1 se of cv
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(GCDlambda[idmin])
  
  idmin <- match(lambda.min, GCDlambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(GCDlambda[idmin])
  
  # return coef for lambda.min
  coef.min = M_coef.powerfamily(newobj, s=lambda.min*nrow(x), type="coefficients")
  
  # return coef for lambda.1se
  coef.1se = M_coef.powerfamily(newobj, s=lambda.1se*nrow(x), type="coefficients")
  
  out <- list(#DrSVM_Fix2.fit = DrSVM.object, 
              lambda = GCDlambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              lambda.min = lambda.min, lambda.1se = lambda.1se,
			  coef.min = t(coef.min), coef.1se = t(coef.1se),
              name = "Misclassification Error"
              )
  class(out) <- "cvs.DrSVM_Fix2"
  out
}
