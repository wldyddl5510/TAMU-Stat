// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// For simplicity, no test data, only training data, and no error calculation.
// X - n x p data matrix
// y - n length vector of classes, from 0 to K-1
// numIter - number of iterations, default 50
// eta - damping parameter, default 0.1
// lambda - ridge parameter, default 1
// beta_init - p x K matrix of starting beta values (always supplied in right format)
// [[Rcpp::export]]
Rcpp::List LRMultiClass_c(const arma::mat& X, const arma::uvec& y, const arma::mat& beta_init,
                          int numIter = 50, double eta = 0.1, double lambda = 1){
    // All input is assumed to be correct
    
    // Initialize some parameters
    int K = max(y) + 1; // number of classes
    int p = X.n_cols;
    int n = X.n_rows;
    arma::mat beta = beta_init; // to store betas and be able to change them if needed
    arma::vec objective(numIter + 1); // to store objective values
    
    // Initialize anything else that you may need
    // indicator matrix
    arma::mat indicator(n, K, arma::fill::zeros);
    indicator.elem((y * n + arma::regspace<arma::uvec>(0, (n - 1)))) = arma::ones<arma::vec>(n);
    
    // taking exp of Xb(n*k mat)
    arma::mat p_mat = arma::exp(X * beta);
    //arma's rowsum: sum(mat, 1)
    arma::colvec rowsum = arma::sum(p_mat, 1);
    // divide row[i](p-vec) by rowsum[i], or equivalently, col[i](n vec) / rowsum (n vec) elementwise
    p_mat.each_col() /= rowsum;
    
    // objective f(beta)
    // log(p_mat)
    arma::mat logp = arma::log(p_mat);
    // add where the indicators are equal to 1
    double logp_sum_nk = arma::accu(logp % indicator);
    //ridge part 
    // square(beta) == (beta % beta) (elementwise)
    double ridge_term = (lambda / 2.0) * arma::accu(square(beta));
    // adding two terms
    objective[0] = (-logp_sum_nk + ridge_term);
    
    // Newton's method cycle - implement the update EXACTLY numIter iterations
    for(int i = 1 ; i <= numIter ; i++) {
        // calculate gradient -> before j loop so that do not have to repeat
        // t(X)(p-indicator) + lambda * beta
        arma::mat gradient = (X.t() * (p_mat - indicator)) + (lambda * beta);
        // arma::mat beta_old = beta;
        for(int j = 0 ; j < K; j++) {
            // index jth column
            arma::colvec p_k_col = p_mat.col(j);
            // get hessian for given column
            // wk, t(X)wX 
            arma::mat loss_hessian = (X.t() * arma::diagmat((p_k_col % (1 - p_k_col))) * X);
            // construct I_p,p
            arma::mat ridge_hessian(p, p, arma::fill::eye);
            // ridge_hessian = lambda * I_p,p
            ridge_hessian *= lambda;
            // sum two mat
            arma::mat hessian_j = (loss_hessian + ridge_hessian);
            // update beta
            beta.col(j) -= (eta * arma::solve(hessian_j, gradient.col(j)));
        }
        // update based on new beta
        // update new p_mat
        // taking exp of Xb(n*k mat)
        p_mat = arma::exp(X * beta);
        //arma's rowsum: sum(mat, 1)
        rowsum = arma::sum(p_mat, 1);
        // divide row[i](p-vec) by rowsum[i], or equivalently, col[i](n vec) / rowsum (n vec) elementwise
        p_mat.each_col() /= rowsum;
        
        // update new objective
        // log(p_mat)
        logp = arma::log(p_mat);
        // main loss
        // add where the indicators are equal to 1
        logp_sum_nk = arma::accu(logp % indicator);
        // ridge part 
        // square(beta) == (beta % beta) (elementwise)
        ridge_term = (lambda / 2.0) * arma::accu(square(beta));
        // adding two terms
        objective[i] = (-logp_sum_nk + ridge_term);
    }
    
    // Create named list with betas and objective values
    return Rcpp::List::create(Rcpp::Named("beta") = beta,
                              Rcpp::Named("objective") = objective);
}