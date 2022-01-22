#' Title
#'
#' @param X mat[n, p]; n row of data with p covariate
#' @param K int; number of clusters
#' @param M mat[K, p]; cluster centers with default NULL; If NULL, initialize K data randomly from X 
#' @param numIter int; number of iteration with default 100
#'
#' @return Y vec[n]; Y[i] is assigned cluster number for ith data
#' @export
#'
#' @examples
#' # Give example
sourceCpp("kmeanscpp.cpp")
MyKmeans_wrapper <- function(X, K, M = NULL, numIter = 100){
  
  n = nrow(X) # number of rows in X
  
  # Check whether M is NULL or not. If NULL, initialize based on K random points from X. If not NULL, check for compatibility with X dimensions.
  X = as.matrix(X)
  # corner case
  if(K == 1) {
    return(rep(1, n))
  }
  p = ncol(X) # p
  # Check whether M is NULL or not. If NULL, initialize based on K randomly selected points from X. If not NULL, check for compatibility with X dimensions.
  # NULL -> init with random k selected X
  # M is NULL
  if(is.null(M)){
    # sample K from n in X's row
    M = X[sample(c(1:n), K), ]
  } else if ((nrow(M) != K) || (ncol(M) != p)) { # compatibility check
    # throw error
    stop("Incompatibility between M and (X or K)")
  }
  # type casting so that we have compatibility
  M = as.matrix(M)
  
  # Call C++ MyKmeans_c function to implement the algorithm
  Y = MyKmeans_c(X, K, M, numIter)
  
  # Return the class assignments
  return(Y)
}