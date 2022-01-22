# Function that implements K-means algorithm. The default number of maximal iterations is 100.
MyKmeans <- function(X, K, M = NULL, numIter = 100) {
  # type casting so that we have compatibility
  X = as.matrix(X)
  X_nrow = nrow(X) # n
  # corner case
  if(K == 1) {
    return(rep(1, X_nrow))
  }
  X_ncol = ncol(X) # p
  # Check whether M is NULL or not. If NULL, initialize based on K randomly selected points from X. If not NULL, check for compatibility with X dimensions.
  # NULL -> init with random k selected X
  # M is NULL
  if(is.null(M)){
    # sample K from n in X's row
    M = X[sample(c(1:X_nrow), K), ]
  } else if ((nrow(M) != K) || (ncol(M) != X_ncol)) { # compatibility check
    # throw error
    stop("Incompatibility between M and (X or K)")
  }
  # type casting so that we have compatibility
  M = as.matrix(M)
  
  # Implement K-means algorithm. It should stop when either (i) the centroids don't change from one iteration to the next, or (ii) the maximal number of iterations was reached, or (iii) one of the clusters has disappeared after one of the iterations (in which case the error message is returned)
  # start iteration
  while(numIter > 0) {
    # get clusters of X's ith rows
    Y = whole_min_cluster(X, M)
    # init new mean mat
    new_mean = matrix(rep(NA, K * X_ncol), K, X_ncol)
    for(k in 1:K) {
      # indices of current cluster
      curr_cluster_indices = which(Y %in% k)
      num_of_curr_cluster = length(curr_cluster_indices)
      # convergence condition -> disappearing cluster
      if(num_of_curr_cluster == 0) {
         stop("a cluster has been disappeared")
      } else if(num_of_curr_cluster == 1) {
        # only one element in the cluster
        new_mean[k, ] = X[curr_cluster_indices, ]
      } else {
        # 2 or more -> vectorization
        new_mean[k, ] = colMeans(X[curr_cluster_indices, ])
      }
    }
    # check convergence condition -> cluster same
    if(identical(M, new_mean)) {
       break;
    }
    # successful updates
    M <- new_mean
    # one iteration
    numIter <- (numIter - 1)
  }
  # Return the vector of assignments Y
  return(Y)
}

# retrieve min dist cluster for given X: (i,j) contains dist(X_i, M_j) -> colmin = colmax(-)
whole_min_cluster <- function(X, M) {
  # (x - y)^2 == x^2 + y^2 - 2xy
  # X^2 is constant in this argmin, so we use X instead to reduce calculations
  mat_dist = outer(rowSums(X), rowSums(M^2), '+') - tcrossprod(X, 2 * M)
  # mimnum col == M_j for argmin given X_i
  return(max.col(-mat_dist))
}