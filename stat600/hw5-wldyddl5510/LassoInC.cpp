#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  // if a > lam -> return a - lam
  if(a > lambda) {
    return (a - lambda);
  }
  // elif a < -lam -> return a + lam
  else if(a < -lambda) {
    return (a + lambda);
  }
  // between [-lam, lam] -> 0
  else {
    return 0;
  }
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  // get num of samples
  int n = Ytilde.n_elem;
  // get residuals
  arma::colvec error_vec = (Ytilde - (Xtilde * beta));
  // take a l2 norm
  // double error_term = std::inner_product(error_vec.begin(), error_vec.end(), error_vec.begin(), 0.0);
  double error_term = arma::dot(error_vec, error_vec) / (2.0 * n);
  // std::cout << error_term;
  // get L1 penalty of beta
  arma::colvec lasso_vec = arma::abs(beta);
  double lasso_term = arma::sum(lasso_vec);
  // return sum
  double total_loss = error_term + (lambda * lasso_term);
  return total_loss;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  int n = Ytilde.n_elem;
  // get p from X
  int p = Xtilde.n_cols;
  // not supplied -> init with zeros
  // if(beta_start == R_NilValue) {
  //  beta_start.zeros();
  // } 
  // get obj value
  double obj_current = lasso_c(Xtilde, Ytilde, beta_start, lambda);
  double obj_gap = obj_current;
  double obj_new;
  arma::colvec beta = beta_start;
  // iterate until conv
  while(obj_gap >= eps) {
    // residual
    arma::colvec r = Ytilde - (Xtilde * beta);
    // new beta
    for(int j = 0 ; j < p; j++) {
      // a is a scalar to take soft
      double a = beta[j] + (arma::dot(Xtilde.col(j), r) / n);
      // soft threshold
      double beta_j_new = soft_c(a, lambda);
      // update residual
      r += ((beta[j] - beta_j_new) * Xtilde.col(j));
      // update beta
      beta[j] = beta_j_new;
    }
    // calculate the new obj
    obj_new = lasso_c(Xtilde, Ytilde, beta, lambda);
    // update obj_gap
    obj_gap = abs(obj_new - obj_current);
    // update obj_value
    obj_current = obj_new;
  }
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  // get n from Ytilde
  // int n = Ytilde.n_elem;
  // get p from Xtilde
  int p = Xtilde.n_cols;

  // get lambda_len
  int n_lambda = lambda_seq.n_elem;
  // init beta_mat: p * n_lambda
  arma::mat beta_mat(p, n_lambda);
  // init value == 0 (p-dim vec)
  arma::colvec init_beta;
  init_beta.zeros(p);
  beta_mat.col(0) = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[0], init_beta, eps);
  // loop if lambda seq is longer than 2
  if(n_lambda >= 2) {
    // loop
    for(int i = 1; i < n_lambda; i++) {
      // warm start
      beta_mat.col(i) = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[i], beta_mat.col(i - 1), eps);
    }
  }
  return beta_mat;
}