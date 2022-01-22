# Use this file to create tests/debug your functions

# Source the functions
source("FunctionsLM.R")

#################################################
# Test "generateY()"
# when X is a vector
# generate beta
p = 1
sigma = 1
beta_1dim = c(rep(3, p))
print(beta_1dim)

# generate X
n = 100
X_1dim_vec = c(rep(2, n))
print(X_1dim_vec)

# generate Y with default seed
Y_1dim_vec = generateY(X_1dim_vec, beta_1dim, sigma)

print(Y_1dim_vec)

# when X is a 1-col matrix
X_1dim_mat = as.matrix(X_1dim_vec)
beta_1dim_mat = as.matrix(beta_1dim)
print(X_1dim_mat)
Y_1dim_mat = generateY(X_1dim_mat, beta_1dim, sigma)
print(Y_1dim_mat)

# when beta is 1-col matrix
Y_1dim_mat_mat = generateY(X_1dim_mat, beta_1dim_mat, sigma)

# since same seed, the result must be same
print(Y_1dim_vec == Y_1dim_mat)
print(Y_1dim_vec == Y_1dim_mat_mat)

# when X is n * p matrix
p = 8 #(n > p)
X_pdim = matrix(rnorm(n * p), n, p)
# random generation of beta
beta_pdim = rnorm(p)
print(beta_pdim)

# generate p_dim data
Y_pdim = generateY(X_pdim, beta_pdim, sigma)
# type testing
Y_pdim_vec = as.vector(Y_pdim)
Y_pdim_mat = as.matrix(Y_pdim)
print(Y_pdim_vec)
print(Y_pdim_mat)

####################################################
# test "calculateBeta()"
# X, Y as a matrix input: X is p dim
beta_LS_mat_mat_pdim = calculateBeta(X_pdim, Y_pdim_mat)
# X mat Y vec
beta_LS_mat_vec_pdim = calculateBeta(X_pdim, Y_pdim_vec)
# X as a vector and Y as a matrix
beta_LS_vec_mat_1dim = calculateBeta(X_1dim_vec, Y_1dim_mat)
# X as a vector and Y as a vector
Y_1dim_vec = as.vector(Y_1dim_vec)
beta_LS_vec_vec_1dim = calculateBeta(X_1dim_vec, Y_1dim_vec)
# X as a matrix and Y as a vector
beta_LS_mat_vec_1dim = calculateBeta(X_1dim_mat, Y_1dim_vec)
# X mat Y mat 1dim
beta_LS_mat_mat_1dim = calculateBeta(X_1dim_mat, Y_1dim_mat)

# compare LS 1, 2, 3, 4
# pdim
print(beta_LS_mat_mat_pdim) # estimate
print(beta_LS_mat_vec_pdim)
print(beta_pdim) # true
# external packages
dimp_model <- lm(Y_pdim_vec ~ X_pdim - 1)
dimp_packages_beta = dimp_model$coef
print(dimp_packages_beta)
# 1dim
print(beta_LS_vec_mat_1dim) # estimates
print(beta_LS_vec_vec_1dim)
print(beta_LS_mat_vec_1dim)
print(beta_LS_mat_mat_1dim)
print(beta_1dim) # true
# external packages
dim1_model <- lm(Y_1dim_vec ~ X_1dim_vec - 1)
dim1_packages_beta = dim1_model$coef
print(dim1_packages_beta)

######################################
# test L2norm ftn
vec1 = c(3,4)
print(L2_norm_square(vec1) == 25)

#####################################
# test calculateEstimateError()
# when input is (scalar, scalar)
beta_scalar = beta_1dim
beta_LS_scalar = as.numeric(beta_LS_vec_vec_1dim)
esti_error_scal_scal = calculateEstimationError(beta_scalar, beta_LS_scalar)

# when imput is (scalar, vector)
beta_LS_1dim_vec = as.vector(beta_LS_vec_vec_1dim)
esti_error_scal_vec = calculateEstimationError(beta_scalar, beta_LS_1dim_vec)

# when input is (scalar, mat)
beta_LS_1dim_mat = as.matrix(beta_LS_vec_vec_1dim)
esti_error_scal_mat = calculateEstimationError(beta_scalar, beta_LS_1dim_mat)

# input: (vec, vec)
beta_vec_1dim = as.vector(beta_1dim)
esti_error_vec_vec_1dim = calculateEstimationError(beta_vec_1dim, beta_LS_1dim_vec)

# input: (mat, mat) 1dim
beta_mat_1dim = as.matrix(beta_1dim)
esti_error_mat_mat_1dim = calculateEstimationError(beta_mat_1dim, beta_LS_1dim_mat)

# when input is (vector, vector)
beta_vec = as.vector(beta_pdim)
beta_LS_vec = as.vector(beta_LS_mat_mat_pdim)
esti_error_vec_vec = calculateEstimationError(beta_vec, beta_LS_vec)

# when input is (vec, mat)
beta_mat = as.matrix(beta_pdim)
esti_error_vec_mat = calculateEstimationError(beta_mat, beta_LS_vec)

# when input is (mat, mat)
beta_LS_mat = as.matrix(beta_LS_mat_mat_pdim)
esti_error_mat_mat = calculateEstimationError(beta_mat, beta_LS_mat)

# 1dim
print(esti_error_scal_scal)
print(esti_error_scal_vec)
print(esti_error_scal_mat)
print(esti_error_vec_vec_1dim)
print(esti_error_mat_mat_1dim)
print(sum(beta_1dim - dim1_packages_beta)^2) # external packages
# pdim
print(esti_error_vec_vec)
print(esti_error_vec_mat)
print(esti_error_mat_mat)
print(sum((beta_pdim - dimp_packages_beta)^2)) # external packages

#################################################
# test calculatePredictionError
# when input is (Y = vec, X = vec, beta_LS = scalar)
prederror_vec_vec_scal = calculatePredictionError(Y_1dim_vec, X_1dim_vec, beta_LS_scalar)

# when input is (Y = vec, X = vec, beta_LS = vec)
prederror_vec_vec_vec = calculatePredictionError(Y_1dim_vec, X_1dim_vec, beta_LS_1dim_vec)

# when input is (Y: vec, X: vec, beta_LS: mat)
prederror_vec_vec_mat = calculatePredictionError(Y_1dim_vec, X_1dim_vec, beta_LS_1dim_mat)

# input: (vec, mat, vec)
prederror_vec_mat_vec = calculatePredictionError(Y_pdim_vec, X_pdim, beta_LS_vec)

# input: (vec, mat, mat)
prederror_vec_mat_mat = calculatePredictionError(Y_pdim_vec, X_pdim, beta_LS_mat)

# input: (mat, mat, vec)
prederror_mat_mat_vec = calculatePredictionError(Y_pdim_mat, X_pdim, beta_LS_vec)

# input: (mat, mat, mat)
prederror_mat_mat_mat = calculatePredictionError(Y_pdim_mat, X_pdim, beta_LS_mat)

# 1dim
print(prederror_vec_vec_scal)
print(prederror_vec_vec_vec)
print(prederror_vec_vec_mat)
print(sum((Y_1dim_mat - fitted(dim1_model))^2)) # external package calculation of pred error
# pdim
print(prederror_vec_mat_vec)
print(prederror_vec_mat_mat)
print(prederror_mat_mat_vec)
print(prederror_mat_mat_mat)
print(sum((Y_pdim - fitted(dimp_model))^2)) # external package calculation of pred error

# first cov -> no exceptions on p == 1?
Xtest_1dim = matrix(rnorm(n * 1), n, 1)
# first cov if Xtest_dim1 == itself
Xtest_first_cov_1dim = Xtest_1dim[, 1, drop = FALSE]
# generate Y for this testset
Ytest_1dim = generateY(Xtest_1dim, beta_1dim_mat, sigma, seed = 678910)
# calculate beta_LS for this test dim1
beta_LS_test_first_cov_1dim = calculateBeta(Xtest_first_cov_1dim, Ytest_1dim)

# dim 10 case
Xtest_pdim = matrix(rnorm(n * 10), n, 10)
# tmp beta
beta_10 = c(rep(5,10))
# first cov of Xtest_dim10
Xtest_first_cov_pdim = Xtest_pdim[, 1, drop = FALSE]
# generate Y for this test dim10
Ytest_pdim = generateY(Xtest_pdim, beta_10, sigma, seed = 678910)
# calculate beta_LS fir this test dim10
beta_LS_test_first_cov_pdim = calculateBeta(Xtest_first_cov_pdim, Ytest_pdim)

# make sure 1dim and pdim shows no different behavior
pred_error_test_first_cov_1dim = calculatePredictionError(Ytest_1dim, Xtest_first_cov_1dim, beta_LS_test_first_cov_1dim)
pred_error_test_first_cov_1dim
pred_error_test_first_cov_pdim = calculatePredictionError(Ytest_pdim, Xtest_first_cov_pdim, beta_LS_test_first_cov_pdim)
pred_error_test_first_cov_pdim
