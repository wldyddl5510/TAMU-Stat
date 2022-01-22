// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//arma::uvec whole_min_cluster_c(const arma::mat& X, const arma::mat& M, int n, int K);

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
        // (x - M)^2 == x^2 + M^2 - 2XM
        // X^2 is constant in this argmin, so we omitted this term to reduce calculations
        // M^2 part
        arma::mat M_square(K, n, arma::fill::zeros);
        // rowSums(M^2)
        arma::colvec M_square_rowsum = arma::sum(arma::square(current_M), 1);
        // copy rowSum to n rows of (rowsum(M^2))
        M_square.each_col() = M_square_rowsum;
        // argmin_{M}(M^2 - 2XM)
        arma::mat mat_dist = M_square.t() - (2 * X) * current_M.t();
        // mimnum col == M_j for argmin given X_i
        Y = arma::index_min(mat_dist, 1);
        
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
        // check convergence condition -> no change in cluster
        // sum(M - N)^2 == 0?
        if(arma::accu(arma::square(current_M - new_mean)) == 0){
            break;
        }
        // successful updates
        current_M = new_mean;
        // one iteration 
        numIter -= 1; // negative one increment
    }
    // add 1 to match index with r
    // Y += 1;
    // Returns the vector of cluster assignments
    return Y;
}