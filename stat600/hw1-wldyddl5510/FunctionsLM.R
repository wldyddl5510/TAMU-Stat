# Generate n-dimensional response Y that follows linear regression model Y = Xbeta + epsilon, where epsilon is normal zero with variance sigma^2 independent across samples. Seed should be set at the beginning of the function
# X - design matrix, rows are n samples
# beta - given parameter vector (could be a vector or a matrix with 1 column)
# sigma - standard deviation of the noise, scalar
# seed  - starting seed value, integer
generateY <- function(X, beta, sigma, seed = 5832652){
  #Set seed and generate Y following linear model
  set.seed(seed)
  # make agnostic to type of X -> convert type to matrix
  X = as.matrix(X)
  # find num of eps (== n)
  n = nrow(X)
  # generate epsilon ~ N(0,sigma)
  eps = rnorm(n, mean = 0, sd = sigma)
  # beta must be agnostic to types -> convert!
  beta = as.matrix(beta)
  # calculate Y
  Y = X %*% beta + eps
  # Return Y
  return(Y)
}

# Calculate beta_LS - least-squares solution, do not use lm function
# X - design matrix, rows are n samples
# Y - response vector (could be a vector or a matrix with 1 column)
calculateBeta <- function(X, Y){
  # Calculate beta_LS
  # make agnostic to input type X and Y
  Y = as.matrix(Y)
  X = as.matrix(X)
  # The solution of linear regression is (X'X)-1X'Y
  X_transpose_X = t(X) %*% X
  projection_mat = solve(X_transpose_X) %*% t(X)
  beta_LS = projection_mat %*% Y
  # Return beta
  return(beta_LS)
}

# Calculate estimation error, defined as ||beta - beta_LS||_2^2
# beta - true coefficient vector (could be a vector or a matrix with 1 column)
# beta_LS - vector estimated by LS (could be a vector or a matrix with 1 column)
calculateEstimationError <- function(beta, beta_LS){
  # Calculate and return error
  # agnostic to type of input
  beta = as.vector(beta)
  beta_LS = as.vector(beta_LS)
  # error as a vector
  error_of_beta = beta - beta_LS
  # L2 norm^2 of error vec
  return(L2_norm_square(error_of_beta))
}


# Calculate prediction error, defined as ||Y - X beta_LS||_2^2
# Y - response vector (could be a vector or a matrix with 1 column)
# X - design matrix, rows are n samples
# beta_LS - vector estimated by LS (could be a vector or a matrix with 1 column)
calculatePredictionError <- function(Y, X, beta_LS){
  # Calculate and return error
  # agostic to type
  beta_LS = as.matrix(beta_LS)
  X = as.matrix(X)
  Y = as.matrix(Y)
  # error by y - yhat
  error_of_yhat = Y - X %*% beta_LS
  # L2 Norm^2 of error vec
  return(L2_norm_square(error_of_yhat))
}

# calculate the L2norm^2 of the given vector
L2_norm_square <- function(x) {
  # L2norm^2 of the vector = square sum of elements
  return(sum(x^2))
}