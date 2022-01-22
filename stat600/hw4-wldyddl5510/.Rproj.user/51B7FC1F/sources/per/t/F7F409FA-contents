###### debugging and testing
# basic dataset
n = 3
p = 4
# X = matrix(c(1:8), nrow = n, ncol = p)
X = matrix(rnorm(n * p), n, p)
X

#### standardization
Xmeans = colMeans(X)
# subtract each col by colmeans
X_centered = X - matrix(Xmeans, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)

# scale
# weights[j] has jth col's norm; (p * 1)
weights = sqrt(colSums(X_centered^2) / n)
# X_centered %*% (1/diag(weights)) in a faster way, by sacrificing the memory
Xtilde = X_centered / rep(weights, rep(n, p))
Xtilde
crossprod(Xtilde)/n

# compare
scale(X_centered)* sqrt(n/(n-1))

#### soft-thesholding
soft1 <- function(a, lambda) {
  return(sign(a) * max(abs(a) - lambda, 0))
}
soft2 <- function(a, lambda) {
  if(a > lambda) {
    return(a - lambda)
  } else if(a < -lambda) {
    return(a + lambda)
  } else {
    return(0)
  }
}

library(microbenchmark)

a = rnorm(1, 2)
lambda = 2.5
soft1(a, lambda) == soft2(a, lambda)
microbenchmark(
  soft1(a, lambda),
  soft2(a, lambda)
)

#### updating
# inside for loop
X = matrix(rnorm(n * p), n , p)
Y = rnorm(n)
beta_t = rnorm(p)
a1 = crossprod(X[ , 1], Y - (X[ , -1] %*% beta_t[-1])) / n
# vectorization
XtXb = crossprod(X) %*% beta_t
Xty = crossprod(X, Y)
first_col_times_excluded_term = XtXb[1] - (sum(X[ , 1]^2) * beta_t[1])
a2 = (Xty[1] - first_col_times_excluded_term) / n
a1
a2

col_calcul_mine <- function(X, Y, beta) {
  XtXb = crossprod(X) %*% beta
  Xty = crossprod(X, Y)
  for(j in 1:ncol(X)) {
    j_col_times_excluded_term = XtXb[j] - (sum(X[ , j]^2) * beta_t[j])
    a = (Xty[j] - j_col_times_excluded_term) / n
  }
}

#### compare speed with partial residual
col_calcul_lec <- function(X, Y, beta) {
  r = Y - X %*% beta
  a = 
  for(j in 1:ncol(X)) {
    a = beta[j] + crossprod(X[ , j], r)
    r = r + X[ , j] * (a - beta[j])
  }
}

library(microbenchmark)

X = matrix(rnorm(400 * 2000), 400, 2000)
Y = rnorm(400)
beta = rnorm(2000)
microbenchmark(
  col_calcul_lec(X, Y, beta),
  col_calcul_mine(X, Y, beta)
)

#### maximum lambda
max(abs(crossprod(X,Y)))


#### return type of list
res = list(c(1,2,3), 10)
A = matrix(rep(0, 3 * 2), 3, 2)
A
res[[2]]
A[ , 1] = res[[1]]
A
b = rep(0, 2)
b
b[1] = res[[2]]
b


#### back scaling
# basic dataset
n = 3
p = 4
# X = matrix(c(1:8), nrow = n, ncol = p)
X = matrix(rnorm(n * p), n, p)
X

#### standardization
Xmeans = colMeans(X)
# subtract each col by colmeans
X_centered = X - matrix(Xmeans, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)

# scale
# weights[j] has jth col's norm; (p * 1)
weights = sqrt(colSums(X_centered^2) / n)

n_lambda = 10

beta_mat_lam = matrix(rep(1, p * n_lambda), p , n_lambda)
beta_mat_lam
weights
beta_mat_lam * weights

#### intercept
Y = rep(1, n)
Y
X = matrix(seq(n * p), n, p)
X
beta_mat_lam = matrix(seq(p * n_lambda), p, n_lambda)
mean(Y) - as.vector(crossprod(colMeans(X), beta_mat_lam))


#### k fold
k = 5
n = 351
table(sample(1:n, size = n) %% k + 1)


#### testing
# parameters
set.seed(1)
n1 = 100
p1 = 2000
n_lambda1 = 30
k1 = 5
eps1 = 0.001

# data
X1 = matrix(rnorm(n1 * p1), n1, p1)
Y1 = rnorm(n1)
# Y2 = rnorm(n1)
# load lib
source("LassoFunctions.R")

# 1. without lambda_seq, fold
res = cvLASSO(X1 ,Y1, lambda_seq = NULL, n_lambda1, k1, fold_ids = NULL, eps1)
res

# 2. properly stop at errors
# 2-1. n mismatch
cvLASSO(rbind(X1, rep(0, p1)), Y1, NULL, n_lambda1, k1, NULL, eps1)

# each function
# 1. standardize
standardized_res = standardizeXY(X1, Y1)
((colSums(standardized_res$Xtilde^2)/n1))
as.numeric(crossprod(standardized_res[[1]])/n1)

sum((colSums((scale(X1) * sqrt(n1/(n1-1)))^2)/n1 - ((colSums(standardized_res$Xtilde^2)/n1)))^2)
standardized_res[[1]]

(scale(Y1, scale = FALSE))[1]
standardized_res[[2]][1]

standardized_res$weights

# 2. soft
soft(1.3, 0.3)
soft(-1.3, 0.3)
soft(0.1, 0.3)
soft(-0.2, 0.3)

# 3. lasso
lasso(standardized_res$Xtilde, standardized_res$Ytilde, rep(1, p1), 1)


# 4. fitLASSOstandardized
ans = fitLASSOstandardized(standardized_res$Xtilde, standardized_res$Ytilde, 0.01, eps = eps1)
ans$fmin
ans$beta

fitLASSOstandardized(standardized_res$Xtilde, standardized_res$Ytilde, 0, eps = 0.00001)

# 5. fitLASSOstandardized_seq
fitlasso_std_res = fitLASSOstandardized_seq(standardized_res$Xtilde, standardized_res$Ytilde, NULL, n_lambda = n_lambda1, eps = eps1)

# 6. fitLASSO
fitlasso_res = fitLASSO(X1, Y1, NULL, n_lambda = n_lambda1, eps = eps1)
all(fitlasso_res$beta_mat == (fitlasso_std_res$beta_mat/standardized_res$weights))
fitlasso_res$beta0_vec
neg_lam_seq = c(-1, 0.02, 3)
fitLASSO(X1, Y1, neg_lam_seq, n_lambda = n_lambda1, eps = eps1)

# plot
colSums(fitlasso_res$beta_mat != 0)
plot(fitlasso_res$lambda_seq, colSums(fitlasso_res$beta_mat != 0))

#### compare
# install.packages("glmnet")
library(glmnet)
scaled_Y = (scale(Y1) * sqrt(n1/(n1-1)))
crossprod(scaled_Y) / n1
X1
scaled_Y
mine = fitLASSO(X1, Y1, NULL, n_lambda = n_lambda1, eps = eps1)
mine_cv = cvLASSO(X1, Y1, NULL, n_lambda = n_lambda1, k = k1, NULL, eps = eps1)
external_res = glmnet(X1, scaled_Y, lambda = mine$lambda_seq, standardize = TRUE, thresh = eps1)
cv_external_res = cv.glmnet(X1, scaled_Y, lambda = mine$lambda_seq, nfolds = 5, foldid = mine_cv$fold_ids)
coef(external_res)
# rbind(mine$beta0_vec, mine$beta_mat)
external_res$beta
mine$beta_mat
mine$beta0_vec

mine_cv$lambda_min
mine_cv$lambda_1se

cv_external_res$lambda.1se
cv_external_res$lambda.min

mine_nonzero = colSums(mine$beta_mat != 0)
plot(mine$lambda_seq, mine_nonzero, col = "red")

external_nonzero = colSums(external_res$beta != 0)
lines(mine$lambda_seq, external_nonzero)

plot(mine_cv$lambda_seq, mine_cv$cvm, col = "red", ylim = c(0.8, 1.5))
lines(mine_cv$lambda_seq, mine_cv$cvm + mine_cv$cvse)
lines(mine_cv$lambda_seq, mine_cv$cvm - mine_cv$cvse)

lines(mine_cv$lambda_seq, cv_external_res$cvm, col = "blue")
lines(mine_cv$lambda_seq, cv_external_res$cvm + mine_cv$cvse, col = "blue")
lines(mine_cv$lambda_seq, cv_external_res$cvm - mine_cv$cvse, col = "blue")
