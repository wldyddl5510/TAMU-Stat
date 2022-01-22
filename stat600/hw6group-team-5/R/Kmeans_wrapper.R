#'  K-means algorithm
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
#' n1 = 40
#' p1 = 20
#' K1 = 5
#' X1 = matrix(rnorm(n1 * p1), n1, p1)
#' MyKmeans(X1, K1, NULL, 100)
#' 
#' n2 = 50
#' p2 = 30
#' K2 = 10
#' X2 = matrix(rnorm(n2 * p2), n2, p2)
#' MyKmeans(X2, K2, NULL, 50)
MyKmeans <- function(X, K, M = NULL, numIter = 100){
  
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
  # adjust index
  Y = Y + 1
  # Return the class assignments
  return(Y)
}