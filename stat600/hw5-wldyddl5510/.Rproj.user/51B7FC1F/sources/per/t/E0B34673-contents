
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################

soft_c(0.5, 0.1)
soft(0.5, 0.1)
soft_c(0.1, 0.5)
soft(0.1, 0.5)
soft_c(-3, 1)
soft(-3, 1)
soft_c(3, 1)
soft(3, 1)

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################

# input 1: large
n1 = 100
p1 = 2000
n_lambda1 = 30
eps1 = 0.001
X1 = matrix(rnorm(n1 * p1), n1, p1)
Y1 = rnorm(n1)
beta_start1 = rep(0, p1)
standardized_res1 = standardizeXY(X1, Y1)
lambda1 = rchisq(1, 3)

# lasso
lasso(standardized_res1$Xtilde, standardized_res1$Ytilde, beta_start1, lambda1)
lasso_c(standardized_res1$Xtilde, standardized_res1$Ytilde, beta_start1, lambda1)

# input 2: small
n2 = 3
p2 = 8
n_lambda2 = 10
eps2 = 0.1
X2 = matrix(seq(n2 * p2), n2, p2)
Y2 = seq(from = n2, to = 1)
beta_start2 = seq(from = p2, to = 1)
standardized_res2 = standardizeXY(X2, Y2)
lambda2 = 0.25

# lasso
lasso(standardized_res2$Xtilde, standardized_res2$Ytilde, beta_start2, lambda2)
lasso_c(standardized_res2$Xtilde, standardized_res2$Ytilde, beta_start2, lambda2)


# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

ans1 = fitLASSOstandardized(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, beta_start1, eps = eps1)
ans1$beta

ans_c1 = fitLASSOstandardized_c(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, beta_start1, eps = eps1)
ans_c1

ans1 = fitLASSOstandardized(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, rep(10, p1), eps = eps1)
ans1$beta

ans_c1 = fitLASSOstandardized_c(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, rep(10, p1), eps = eps1)
ans_c1

ans2 = fitLASSOstandardized(standardized_res2$Xtilde, standardized_res2$Ytilde, lambda2, beta_start2, eps = eps2)
ans2$beta

ans_c2 = fitLASSOstandardized_c(standardized_res2$Xtilde, standardized_res2$Ytilde, lambda2, beta_start2, eps = eps2)
ans_c2

ans2 = fitLASSOstandardized(standardized_res2$Xtilde, standardized_res2$Ytilde, lambda1, rep(10, p2), eps = eps2)
ans2$beta

ans_c2 = fitLASSOstandardized_c(standardized_res2$Xtilde, standardized_res2$Ytilde, lambda1, rep(10, p2), eps = eps2)
ans_c2

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

library(microbenchmark)

microbenchmark(
  fitLASSOstandardized(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, beta_start1, eps = eps1),
  fitLASSOstandardized_c(standardized_res1$Xtilde, standardized_res1$Ytilde, lambda1, beta_start1, eps = eps1)
)

# Unit: milliseconds
# expr
# fitLASSOstandardized(standardized_res1$Xtilde, standardized_res1$Ytilde,      lambda1, beta_start1, eps = eps1)
# fitLASSOstandardized_c(standardized_res1$Xtilde, standardized_res1$Ytilde,      lambda1, beta_start1, eps = eps1)
# min       lq     mean   median       uq      max neval
# 23.0175 36.91645 51.27286 44.60145 65.66315 127.6130   100
# 1.1148  1.87370  2.57820  2.51125  3.10585   5.5865   100

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

seq_ans1 = fitLASSOstandardized_seq(standardized_res1$Xtilde, standardized_res1$Ytilde, NULL, n_lambda = n_lambda1, eps = eps1)
mine_seq_ans1 = fitLASSOstandardized_seq_c(standardized_res1$Xtilde, standardized_res1$Ytilde, seq_ans1$lambda_seq, eps = eps1)
seq_ans1$beta_mat == mine_seq_ans1

seq_ans2 = fitLASSOstandardized_seq(standardized_res2$Xtilde, standardized_res2$Ytilde, NULL, n_lambda = n_lambda1, eps = eps2)
mine_seq_ans2 = fitLASSOstandardized_seq_c(standardized_res2$Xtilde, standardized_res2$Ytilde, seq_ans2$lambda_seq, eps = eps2)
seq_ans2$beta_mat == mine_seq_ans2

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

microbenchmark(
  fitLASSOstandardized_seq(standardized_res1$Xtilde, standardized_res1$Ytilde, seq_ans1$lambda_seq, n_lambda = n_lambda1, eps = eps1),
  fitLASSOstandardized_seq_c(standardized_res1$Xtilde, standardized_res1$Ytilde, seq_ans1$lambda_seq, eps = eps1)
)

# Unit: milliseconds
# expr
# fitLASSOstandardized_seq(standardized_res1$Xtilde, standardized_res1$Ytilde,      seq_ans1$lambda_seq, n_lambda = n_lambda1, eps = eps1)
# fitLASSOstandardized_seq_c(standardized_res1$Xtilde, standardized_res1$Ytilde,      seq_ans1$lambda_seq, eps = eps1)
# min        lq       mean    median        uq       max neval
# 1249.4738 1556.3698 2191.61770 2016.8261 2676.9973 3785.3308   100
# 41.1823   64.9261   94.67703   97.8846  121.4509  173.6647   100

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)
out2 = fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq)
# compare result
which(outl$beta != out2)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)

# Unit: milliseconds
# expr       min        lq
# fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq) 2937.3759 3177.9354
# fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq)   67.4634   75.3767
# mean     median        uq       max neval
# 4035.69658 3547.24685 4950.4095 5906.8534    10
# 92.24669   93.94165   99.1957  130.6484    10
