# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # get the length param
  n = length(Y)
  p = ncol(X)
  # [ToDo] Center Y
  Ymean = mean(Y)
  Ytilde = Y - Ymean
  # [ToDo] Center and scale X
  # center X : (p * 1)
  Xmeans = colMeans(X)
  # subtract each col by colmeans
  X_centered = X - matrix(Xmeans, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  
  # scale
  # weights[j] has jth col's norm; (p * 1)
  weights = sqrt(colSums(X_centered^2) / n)
  # X_centered %*% (1/diag(weights)) in a faster way, by sacrificing the memory
  Xtilde = X_centered / rep(weights, rep(n, p))
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  # if-else statement seems faster than one line function on lecture notes
  # if-else statement would not work on vectors, but only testing scalar
  # large positive > lambda
  if(a > lambda) {
    # sgn(a) * (abs(a) - lambda) = a - lambda
    return(a - lambda)
  } else if(a < -lambda) { # smaller negative < lambda
    # sgn(a) * (abs(a) - lambda) = -(-a - lambda) = a + lambda
    return(a + lambda)
  } else { # -lambda <= a <= lambda
    return(0)
  }
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  # ||Y - Xb||^2 = crossprod(Y - Xb)
  # get num of sample
  n = length(Ytilde)
  # L2 error estimating term
  error_term = sum((Ytilde - (Xtilde %*% beta))^2) / (2 * n)
  # L1 penalty of beta
  lasso_term = sum(abs(beta))
  # total loss = sum with lambda on lasso term
  total_loss = error_term + (lambda * lasso_term)
  return(total_loss)
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  # get n from Y
  n = length(Ytilde)
  # check compatibility with X row
  if(nrow(Xtilde) != n) {
    stop("number of data in X does not match with number of data in Y")
  }
  #[ToDo]  Check that lambda is non-negative
  if(lambda < 0) {
    stop("Lambda must be non-negative")
  }
  #[ToDo]  Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  # get p from X
  p = ncol(Xtilde)
  # Non supplied -> init with zeros
  if(is.null(beta_start)) {
    beta_start = rep(0, p)
  } else{ # check compatibility of beta
    if(length(beta_start) != p) {
      stop("p(= length of beta) must be equal to column of X")
    }
  }
  
  #[ToDo]  Coordinate-descent implementation. Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999, your should return fmin = 0.99999, and not have another iteration
  
  # get init obj value
  obj_current = lasso(Xtilde, Ytilde, beta_start, lambda)
  # obj_gap is abs(obj[t] - obj[t-1])
  # init obj_gap is init value of obj function in beta_start
  obj_gap = obj_current
  # break when obj_gap < eps
  # store beta
  beta = beta_start
  while(obj_gap >= eps) {
    # store repeating computations: partial residual
    r = Ytilde - (Xtilde %*% beta)
    # Update beta's
    for(j in 1:p) {
      # a is a vector to take soft
      a = as.numeric(beta[j] + (crossprod(Xtilde[ , j], r) / n))
      # soft threshold
      beta_j_new = soft(a, lambda)
      # update residual
      r = r + ((beta[j] - beta_j_new) * Xtilde[ , j])
      # update beta
      beta[j] = beta_j_new
    }
    # calculate the new obj
    obj_new = lasso(Xtilde, Ytilde, beta, lambda)
    # update obj gap
    obj_gap = abs(obj_new - obj_current)
    # update obj value
    obj_current = obj_new
  }
  # last updated obj value = fmin
  fmin = obj_current
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  # get n from Y and p from X
  n = length(Ytilde)
  p = ncol(Xtilde)
  # check compatibility with X row
  if(nrow(Xtilde) != n) {
    stop("number of data in X does not match with number of data in Y")
  }
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0, and make sure the values are sorted from largest to smallest. If none of the supplied values satisfy the requirement, print the warning message and proceed as if the values were not supplied.
  # lambda_seq is not null
  if(!is.null(lambda_seq)) {
    # indexing values >= 0
    # explicit_num_lambda = length(lambda_seq)
    lambda_seq = lambda_seq[which(lambda_seq >= 0)]
    # check num of non-neg
    if(length(lambda_seq) == 0) {
      warning("no non-neg lambda's. lambda_seq will be regarded as NULL")
      # lambda_max >= max_j abs(t(Xj)Y) / n
      lambda_max = max(abs(crossprod(Xtilde, Ytilde)))/n
      lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
    }
  } else { # lambda seq is null
    # If lambda_seq is not supplied, calculate lambda_max (the minimal value of lambda that gives zero solution), and create a sequence of length n_lambda as
    lambda_max = max(abs(crossprod(Xtilde, Ytilde)))/n
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  # sort
  lambda_seq = sort(lambda_seq, decreasing = TRUE)
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda (make sure supplied eps is carried over). Use warm starts strategy discussed in class for setting the starting values.
  # get true length
  true_lambda_len = length(lambda_seq)
  # init beta and fmin
  beta_mat = matrix(0, p, true_lambda_len)
  fmin_vec = rep(0, true_lambda_len)
  
  # perform init step
  beta_fmin_list = fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[1], beta_start = beta_mat[ , 1], eps)
  beta_mat[ , 1] = beta_fmin_list$beta
  fmin_vec[1] = beta_fmin_list$fmin
  # loop from 2 to lambda_len if true_lambda_len >= 2
  if(true_lambda_len >= 2) {
    for(i in 2:true_lambda_len) {
      # warm start -> start from beta_mat[ , i - 1]: res of prev lambda
      beta_fmin_list = fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_start = beta_mat[ , (i - 1)], eps)
      beta_mat[ , i] = beta_fmin_list$beta
      fmin_vec[i] = beta_fmin_list$fmin
    }
  }
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardizing = standardizeXY(X, Y)
  Xtilde = standardizing$Xtilde
  Ytilde = standardizing$Ytilde
  Ymean = standardizing$Ymean
  Xmeans = standardizing$Xmeans
  weights = standardizing$weights
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
  lam_beta_fmin_list = fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector for each lambda
  # back scaling -> just multiply weights, columnwise
  # beta_mat[j] = beta[j](p*1) / weights(p*1) elementwise
  beta_mat = lam_beta_fmin_list$beta_mat / weights
  
  # get beta0
  beta0_vec = Ymean - as.vector(crossprod(Xmeans, beta_mat))
  
  # store true lambda
  lambda_seq = lam_beta_fmin_list$lambda_seq
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # get num of samples
  n = length(Y)
  # [ToDo] Fit Lasso on original data using fitLASSO
  lam_beta_beta0_list = fitLASSO(X, Y, lambda_seq, n_lambda, eps)
  
  # get actual num of lambda and lambda (non-neg)
  lambda_seq = lam_beta_beta0_list$lambda_seq
  actual_num_of_lam = length(lambda_seq)
  
  # save other outputs
  beta_mat = lam_beta_beta0_list$beta_mat
  beta0_vec= lam_beta_beta0_list$beta0_vec
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds. If fold_ids is not NULL, split the data according to supplied fold_ids.
  if(is.null(fold_ids)) {
    # make sure there exists a positive k
    if(k <= 0) {
      stop("either fold_ids or k > 0 must be provided.")
    }
    fold_ids = (sample(1:n, size = n) %% k) + 1
  } else {
    # get the number of fold from provided fold
    k = max(folds_id)
  }
  # X[fold_ids, ] and Y[fold_ids] is splited data
  
  # cvm and cvse
  # cvm = rep(0, actual_num_of_lam)
  # cvse = rep(0, actual_num_of_lam)
  # cv for each fold
  cv_fold = matrix(0, k, actual_num_of_lam)
  # [ToDo] Calculate LASSO on each fold using fitLASSO, and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  for(fold in 1:k) {
    # create training and validation set
    X_train = X[fold_ids != fold, ]
    Y_train = Y[fold_ids != fold]
    X_val = X[fold_ids == fold, ]
    Y_val = Y[fold_ids == fold]
    
    # fit lasso on training set
    fold_lam_beta_beta0_list = fitLASSO(X_train, Y_train, lambda_seq, eps)
    
    # cv for each fold in validation set
    # get length of validation set
    num_val = length(Y_val)
    # intercept = (num_val * n_lambda) vec
    intercept_by_row = rep(fold_lam_beta_beta0_list$beta0_vec, rep(num_val, actual_num_of_lam))
    # cv each fold (1 * n_lambda) length
    cv_fold[fold, ] = colMeans((Y_val - intercept_by_row - (X_val %*% fold_lam_beta_beta0_list$beta_mat))^2)
  }
  # get cvm and cvse
  # cvm = vector length = n_lambda
  cvm = colMeans(cv_fold)
  # cvse = (1 * n_lambda)
  cvse = apply(cv_fold, 2, sd) / sqrt(k)
  
  # [ToDo] Find lambda_min
  # get min lambda -> if multiple exists only first(=largest) lambda
  min_cv_lam = which.min(cvm)
  min_cv = cvm[min_cv_lam]
  lambda_min = lambda_seq[min_cv_lam]
  # [ToDo] Find lambda_1SE
  cv_se_min_lam = cvse[min_cv_lam]
  cv_criterion = min_cv + cv_se_min_lam
  # since lambda_seq is sorted, get the first element
  largest_lam_given_criterion = which(cvm <= cv_criterion)[1]
  lambda_1se = lambda_seq[largest_lam_given_criterion]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}