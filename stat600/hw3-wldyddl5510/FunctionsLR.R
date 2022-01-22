# Function that implements multi-class logistic regression.
#############################################################
# Description of supplied parameters:
# X - n x p training data, 1st column should be 1s to account for intercept
# y - a vector of size n of class labels, from 0 to K-1
# Xt - ntest x p testing data, 1st column should be 1s to account for intercept
# yt - a vector of size ntest of test class labels, from 0 to K-1
# numIter - number of FIXED iterations of the algorithm, default value is 50
# eta - learning rate, default value is 0.1
# lambda - ridge parameter, default value is 0.1
# beta_init - (optional) initial starting values of beta for the algorithm, should be p x K matrix 

## Return output
##########################################################################
# beta - p x K matrix of estimated beta values after numIter iterations
# error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
# error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
# objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
LRMultiClass <- function(X, y, Xt, yt, numIter = 50, eta = 0.1, lambda = 1, beta_init = NULL){
  ## Check the supplied parameters as described. You can assume that X, Xt are matrices; y, yt are vectors; and numIter, eta, lambda are scalars. You can assume that beta_init is either NULL (default) or a matrix.
  ###################################
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  if(any(X[ , 1] != 1)) {
    stop("1st column of X must be 1")
  }
  if(any(Xt[ , 1] != 1)) {
    stop("1st column of Xt must be 1")
  }
  # Check for compatibility of dimensions between X and Y
  # num of sample
  n = nrow(X)
  # match with y
  if(n != length(y)) {
    stop("Num of samples in X and Y mismatch")
  }
  # Check for compatibility of dimensions between Xt and Yt
  # num of test samples
  ntest = nrow(Xt)
  # match with Yt
  if(ntest != length(yt)) {
    stop("Num of samples in Xtest and Ytest mismatch")
  }
  # Check for compatibility of dimensions between X and Xt
  # num of covariates
  p = ncol(X)
  if(p != ncol(Xt)) {
    stop("Num of covariates in X and Xt mismatch")
  }
  # Check eta is positive
  if(eta <= 0) {
    stop("Eta must be positive")
  }
  # Check lambda is non-negative
  if(lambda < 0) {
    stop("Lambda must be nonnegative")
  }
  # Check whether beta_init is NULL. If NULL, initialize beta with p x K matrix of zeroes. If not NULL, check for compatibility of dimensions with what has been already supplied.
  # k = maximum class
  k = max(y) + 1
  if(is.null(beta_init)) {
    # init beta
    beta_init = matrix(0, p, k)
  } else {
    if ((nrow(beta_init) != p) || (ncol(beta_init) != k)) {
      stop("dim of beta not compatible")
    }
  }
  ## Calculate corresponding pk, objective value f(beta_init), training error and testing error given the starting point beta_init
  ##########################################################################
  error_train = rep(0, numIter + 1)
  error_test = rep(0, numIter + 1)
  objective = rep(0, numIter + 1)
  
  # adjust y to be [1,k]
  y = y + 1
  yt = yt + 1
  # calculate init value -> training set
  beta = beta_init
  # p_ij = ith sample's prob of in jth class
  p_mat_training = p_mat_generator(X, beta)
  # objective f(beta)
  objective[1] = multi_logit_ridge_objective(beta, p_mat_training, X, y, lambda)
  # estimation
  # calculate the error == ratio of mismatch
  error_train[1] = error_calculator(y, p_mat_training, n)
  
  # test set
  p_mat_test = p_mat_generator(Xt, beta)
  error_test[1] = error_calculator(yt, p_mat_test, ntest)
  
  ## Newton's method cycle - implement the update EXACTLY numIter iterations
  ##########################################################################
  for(i in 2:(numIter + 1)) {
    
    # Within one iteration: perform the update, calculate updated objective function and training/testing errors in % 
    # update beta
    # iterate k and update beta
    
    gradient = multinom_logit_gradient(X, y, p_mat_training, beta, n, k, lambda)
    beta_old = beta
    for(j in 1:k) {
      # kth column of p_mat
      p_k_col = as.vector(p_mat_training[ , j])
      hessian = multinom_logit_hessian_col(X, p_k_col, p, lambda)
      # update beta
      beta[, j] = beta[, j] - (eta * solve(hessian, gradient[, j]))
    }
    # update based on new beta
    p_mat_training = p_mat_generator(X, beta)
    p_mat_test = p_mat_generator(Xt, beta)
    objective[i] = multi_logit_ridge_objective(beta, p_mat_training, X, y, lambda)
    error_train[i] = error_calculator(y, p_mat_training, n)
    error_test[i] = error_calculator(yt, p_mat_test, ntest)
  }
  
  
  ## Return output
  ##########################################################################
  # beta - p x K matrix of estimated beta values after numIter iterations
  # error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
  # error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
  # objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
  return(list(beta = beta, error_train = error_train, error_test = error_test, objective =  objective))
}

# objective function
# @params: p: p_matrix(n,p); X: X matrix(n,p); y_adjusted: result class vector(n) -> adjusted by [1, k]; 
# @params: lambda: scalar
# @returns: scalar
multi_logit_ridge_objective <- function(beta, p_mat, X, y_adjusted, lambda) {
  # \sum_k 1_{y=k} log p as a matrix form
  logp = log(p_mat)
  # t(X)[seq(0, nrow(X) - 1) * ncol(X) + (y + 1)] is (n * 1) vector indexing p_k(i) which is y = k
  logp_sum_nk = sum(t(logp)[seq(0, nrow(logp) - 1) * ncol(logp) + y_adjusted])
  # ridge part
  ridge_term = (lambda / 2) * sum(beta^2)
  # adding two terms
  return(-logp_sum_nk + ridge_term)
}

# calculating gradient
# @params: X: data matrix(n,p); y_adjusted: class vector(n) adjusted [1, k], p_mat: n*k prob matrix
# @returns: p*k matrix: each column has kth class'es gradient
multinom_logit_gradient <- function(X, y_adjusted, p_mat, beta, n, k, lambda) {
  # creating (1_{Y = l})_{l = 1 to k}: n * k mat
  indi = matrix(0, nrow = n, ncol = k)
  indi[cbind(seq_along(y_adjusted), y_adjusted)] <- 1
  # (t(X)(p - indi)) + lambda * beta
  return(crossprod(X, p_mat - indi) + (lambda * beta))
}

# columnwise -> for loop version
multinom_logit_hessian_col <- function(X, p_k_col, p, lambda) {
  # wk, t(X)wX
  loss_hessian = crossprod(X, (p_k_col * (1-p_k_col)) * X)
  # lambda I: p * p
  ridge_hessian = lambda * diag(nrow = p)
  # adding two
  return(loss_hessian + ridge_hessian)
}

# p_mat from X and beta
# @params: X: data matrix (n*p); beta: coeff matrix (p*k)
# @returns: n * k matrix prob of i row sample in k col class
p_mat_generator <- function(X, beta) {
  # X^t b : (n * k) mat
  Xb = X %*% beta
  # taking exponetional on (n * k) mat is also (n * k) matrix
  exp_Xb = exp(Xb)
  # summing up by row yields (n * 1): each as a all prob for the sample i -> denominator
  # p_ij = ith sample's prob of in jth class
  return(exp_Xb / rowSums(exp_Xb))
}

# calculating error using y and p
# derive yhat from p and calculate y != yhat
# @params: y_adjusted: n vector from [1, k] (adjusted) ; p_mat: (n * k) prob matrix
# @returns: error (percentized)
error_calculator <- function(y_adjusted, p_mat, n) {
  # max prob element is yhat -> -1 to get [0,k-1] from [1,k]
  y_hat = max.col(p_mat, ties.method = "first")
  #print(y_hat)
  #print(y_adjusted)
  # count num of mismatch
  # * 100 to percentize
  return((sum(y_hat != y_adjusted) / n) * 100)
}
