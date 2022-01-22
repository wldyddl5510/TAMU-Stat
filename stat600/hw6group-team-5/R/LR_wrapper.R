
#' Multi logistic regression with ridge penalty 
#'
#' @param X mat[n, p]: n row of data, p col of covariates
#' @param y vec[n]: each element is lie on 0 ~ K-1, denoting the class of the sample
#' @param numIter number of iterations, default 50
#' @param eta damping parameter, default 0.1, must be positive
#' @param lambda ridge parameter, default 1, must be non negative
#' @param beta_init mat[p, k]: starting beta values, if null, initialize with 0 matrix.
#'
#' @return out - list of updated beta and objective after iterations
#' @export
#'
#' @examples 
#' n1 = 50
#' p1 = 30
#' k1 = 10
#' X1 = matrix(rnorm(n1, (p1 - 1)), n1, p1 - 1)
#' X1 = cbind(rep(1, n1), X1)
#' y1 = sample(k1, n1, replace = TRUE) - 1
#' LRMultiClass(X1, y1, 50, 0.1, 1, NULL)
#' 
#' beta_init = matrix(1, p1, k1)
#' LRMultiClass(X1, y1, 100, 0.5, 0.5, beta_init)
LRMultiClass <- function(X, y, numIter = 50, eta = 0.1, lambda = 1, beta_init = NULL){
  
  # Compatibility checks from HW3 and initialization of beta_init
  # convert to matrix so that no compatibility miss
  X = as.matrix(X)
  # Check that the first column of X is 1's, if not - display appropriate message and stop execution.
  if(any(X[ , 1] != 1)) {
    stop("1st column of X must be 1")
  }
  # Check for compatibility of dimensions between X and Y
  # num of sample
  n = nrow(X)
  # match with y
  if(n != length(y)) {
    stop("Num of samples in X and Y mismatch")
  }
  # num of covariates
  p = ncol(X)
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
  
  # Call C++ LRMultiClass_c function to implement the algorithm
  out = LRMultiClass_c(X, y, beta_init, numIter, eta, lambda)
  
  # Return the class assignments
  return(out)
}
