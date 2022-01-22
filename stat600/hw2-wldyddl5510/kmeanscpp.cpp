// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

arma::uvec whole_min_cluster_c(const arma::mat& X, const arma::mat& M, int n, int K);

// X - n x p matrix
// K - number of clusters
// M - K x p cluster centers (always given)
// numIter - maximal number of iterations
// [[Rcpp::export]]
arma::uvec MyKmeans_c(const arma::mat& X, int K, const arma::mat& M, int numIter = 100){
  // All input is assumed to be correct
  // Initialize some parameters
  int n = X.n_rows;
  int p = X.n_cols;
  arma::uvec Y(n); // to store cluster assignments
  
  // Initialize any additional parameters if needed
  arma::mat current_M = M;
  
  // For loop with kmeans algorithm
  while(numIter > 0) { 
    // get clusters of X's ith rows
    Y = whole_min_cluster_c(X, current_M, n, K);
    // init new mean mat
    arma::mat new_mean(K, p);
    
    // iterate each cluster
    for(int k = 0 ; k < K ; k++) {
      // indices of current cluster
      arma::uvec curr_cluster_indices = arma::find(Y == k);
      // get number of elements in current cluster
      int num_of_curr_cluster = curr_cluster_indices.n_elem;
      // convergence condition -> disappearing cluster
      if(num_of_curr_cluster == 0) {
        Rcpp::stop("a cluster has been disappeared");
      } else if(num_of_curr_cluster == 1) {
        // only one element in the cluster
        new_mean.row(k) = X.rows(curr_cluster_indices);
      } else {
        // 2 or more -> new mean = rowMeans(current cluster data)
        new_mean.row(k) = arma::mean(X.rows(curr_cluster_indices), 0);
      }
    }
    // check convergence condition -> cluster same
    // sum(M - N)^2 == 0?
    if(arma::accu(arma::square(current_M - new_mean)) == 0){
      break;
    }
    // successful updates
    current_M = new_mean;
    // one iteration
    numIter -= 1;
  }
  
  // Returns the vector of cluster assignments
  return Y;
}

// retrieve min dist cluster for given X: (i,j) contains dist(X_i, M_j)
/*
 * @params: mat X[n, p]; mat M[K, p] - ith row is ith cluster's center; int n - num of data ; int k - num of clusters
 * @returns: uvec y[n]: ith element is the number of cluster for ith data.
 */
arma::uvec whole_min_cluster_c(const arma::mat& X, const arma::mat& M, int n, int K) {
  // (x - y)^2 == x^2 + y^2 - 2xy
  // X^2 is constant in this argmin, so we use X instead to reduce calculations
  // M^2 part
  arma::mat M_square(K, n);
  // rowSums(M^2)
  arma::colvec M_square_rowsum = arma::sum(arma::square(M), 1);
  // copy rowSum to n rows of (rowsum(M^2))
  M_square.each_col() = M_square_rowsum;
  // argmin_{m}(m^2 - 2xm)
  arma::mat mat_dist = M_square.t() - (2 * X) * M.t();
  // mimnum col == M_j for argmin given X_i
  return arma::index_min(mat_dist, 1);
}

// [[Rcpp::export]]
arma::rowvec test_rows(arma::mat& X, arma::uvec& indices) {
  int num_of_curr_cluster = indices.n_elem;
  // convergence condition -> disappearing cluster
  if(num_of_curr_cluster == 0) {
    Rcpp::stop("a cluster has been disappeared");
  } else if(num_of_curr_cluster == 1) {
    // only one element in the cluster
    return X.rows(indices);
  } else {
    // 2 or more -> new mean = rowMeans(current cluster data)
    return arma::mean(X.rows(indices), 0);
  }
}