cvs.DrSVM_Fix21 <- function(x, y, lambda2, lambda = NULL, nfolds = 5, foldid, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  y <- c(-1, 1)[as.factor(drop(y))]
  
  n <- nrow(x)
  p <- ncol(x)
  
  source("D:\\GitHub\\powerfamily\\DrSVM\\M_DrSVM_Fix2.r")
  source("D:\\GitHub\\powerfamily\\O_utilities.r")
  source("D:\\GitHub\\powerfamily\\M_coef.GCDpower.r")
  source("D:\\GitHub\\powerfamily\\M_p.GCDpower.r")
  
  x <- as.matrix(x)
  if(missing(lambda2))
  {
    lambda2 = 1 * n
  }
  
  if(!is.null(lambda))
  {
    GCDlambda = lambda
  } else
  {
  # Fit a model first, in order to get lambda sequence
    GCDpower.object <- GCDpower(x, y, delta = 0.01, qv = 2, standardize=F)
    GCDlambda <- GCDpower.object$lambda
  }
  
  DrSVM.object <- DrSVM_Fix2(x, y, lambda2 = lambda2 * n, ...)
  
  lambda <- DrSVM.object$lambda
  
  # record the number of non-zero coefficients for each lambda
  nbeta <- t(DrSVM.object$beta)
  
  dimnames(nbeta) <- list(NULL, NULL)
  nzel <- function(x, which) if (any(x)) 
    which[x] else NULL
  t.beta <- abs(as.matrix(t(nbeta))) > 0
  
  nz <- (if (ncol(nbeta) == 1) 
    apply(t.beta, 2, nzel, seq(p)) else apply(t.beta, 1, nzel, seq(p)))
  nz <- sapply(nz, length)
  names(nz) <- paste("s", seq(nz), sep = "")
  
  # Start to initialize folder
  # set.seed(124)
  foldid <- sample(rep(seq(nfolds), length = n)) 
  outlist <- as.list(seq(nfolds))
  
  # Within each folder, fit a model using training set
  for (i in seq(nfolds)) {
    which <- foldid == i
	x_sub = x[!which, , drop = FALSE]
    outlist[[i]] <- DrSVM_Fix2(x = x_sub, 
                               y = y[!which], 
                               lambda2 = lambda2 * nrow(x_sub), ...)
  }
  
  # Compute the misclassification rate
  predmat <- matrix(NA, length(y), length(GCDlambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    obj <- NULL
    obj$b0 <- (apply(fitobj$beta0, 1, mean))
    obj$beta <- t(fitobj$beta)
    obj$lambda <- fitobj$lambda1
    preds <- predict.GCDpower(obj, x[which, , drop = FALSE],
                              s=GCDlambda*nrow(x_sub), type="class")
    predmat[which, seq(GCDlambda)] <- preds
    nlams[i] <- length(GCDlambda)
  }
  # nfold rows and l column
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
  DrSVM.object = newobj
  
  # For each folder, record the length of lambda
  N <- apply(good, 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  out <- list(lambda = GCDlambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              nzero = nz, name = "Misclassification Error"
              ,DrSVM_Fix2.fit = DrSVM.object
  )
  #lamin <- getmin(lambda, cvm, cvsd)
  
  # Get the appropriate lambda
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(GCDlambda[idmin])
  idmin <- match(lambda.min, GCDlambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(GCDlambda[idmin])
  lamin <- list(lambda.min = lambda.min, lambda.1se = lambda.1se)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.DrSVM_Fix2"
  obj
}

