ridge_fit_ <- function(d, targ_col, Embedding,theta, lambda,alp){
  #####################################################################
  ###### Explanation of the algorithm to recover the S-map coefficients
  Edim <- length(Embedding)
  if (ncol(d) > Edim){
    d <- d[,Embedding]
  }
  coeff_names <- sapply(colnames(d),function(x) paste("d", targ_col, "d", x, sep = ""))
  #### This is a Nx(E +1) matrix where the first column is the target column
  block <- cbind(d[2:dim(d)[1],targ_col],d[1:(dim(d)[1]-1),])
  #### This is a (E + 1) vector with the sd of the time series
  norm_consts <- apply(block, 2, function(x) sd(x))
  #### This is the rescaled block matrix so that the time series are centered around 0
  block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))
  
  ##### Sequence from 1 to N
  lib <- 1:dim(block)[1]
  pred <- 1:dim(block)[1]
  
  #### NxE matrix of zeros
  coeff <- array(0,dim=c(length(pred),Edim + 1))
  colnames(coeff) <- c('c0', coeff_names)
  #### These are going to be the elements of the Jacobian
  coeff <- as.data.frame(coeff)
  
  
  
  lm_regularized <- function(y, x, ws, lambda, dimension, subset = seq_along(y)){
    #### y is the target colum (N-1) vector.
    #### x are all the others (N-1) x E vector
    #### Ws is a (N-1) vector
    x <- x[subset,]
    y <- y[subset]
    ws <- ws[subset]
    #########################################################
    WWs = diag(ws)
    Xx = as.matrix(x)
    #### For constant term in linear fit
    Xx = cbind(1, Xx)
    coeff <- solve(t(Xx) %*% WWs %*% Xx + lambda*nrow(Xx)*diag(1,dimension + 1)) %*% t(Xx) %*%(ws * y)
    coeff <- t(coeff)
    
    return(coeff)
  }
  
  ############### Remember: block[i,1] is the i-th row of the target column 
  ###############           block[i,>1] is the i-th row of any other column
  ############### pred and lib are sequences from 1 to N
  fit_error = rep(0,length(pred))
  for (ipred in 1:length(pred)){
    #target point is excluded from the fitting procedure
    libs = lib[-pred[ipred]]
    # q is a (N-1)xE matrix with the ipred-th entry of the rescaled time series (yes it is a N-1 repetation of one value)
    q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
                ncol=Edim, nrow=length(libs), byrow = T)
    ##########################################################
    ######### Here compute the weigths. Wx is going to be a (N-1) vector of weights
    distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
    Krnl = match.fun(Regression.Kernel)
    Ws = Krnl(distances, theta)
    ##########################################################
    ####### svd_fit gives me the Jacobian element at each time step ipred
    fit <- lm_regularized(block[libs,1],block[libs,2:dim(block)[2]],Ws, lambda, Edim)
    coeff[ipred,] <- fit
  }
  return(coeff)
  
}

ridge_fit <- cmpfun(ridge_fit_)