newobj <- function(DrSVMobject)
{
  index = 0
  if(sum(diff(DrSVMobject$lambda1) == 0) > 0)
  {
    index = which(diff(DrSVMobject$lambda1) == 0)
  }
  object = NULL
  object$b0 <- apply(DrSVMobject$beta0, 1, mean)
  names(object$b0) <- paste("s", seq(object$b0), sep="")
  
  object$beta <- t(DrSVMobject$beta)
  dimnames(object$beta) <- list(paste("V", seq(nrow(object$beta)), sep=""),
                                paste("s", (seq(ncol(object$beta))-1), sep=""))
  object$lambda <- DrSVMobject$lambda1
  if(sum(index) > 0)
  {
    object$b0 <- object$b0[-index]
	object$beta <- object$beta[, -index]
	object$lambda <- object$lambda[-index]
  }
  
  object
}

# This code is correctly plugging in coefficients for lambda!
M_coef.powerfamily <- function(object, s = NULL, type = c("coefficients", 
                                                          "nonzero"), ...) {
  type <- match.arg(type)
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    # Plug-In
	if(length(s) == 1)
    {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
        nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else
    {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] %*% diag(lamlist$frac) + 
        nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    
    
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  if (type == "coefficients") 
    return(nbeta)
  if (type == "nonzero") 
    return(nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE))
} 