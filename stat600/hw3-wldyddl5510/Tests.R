# This is a script to save your own tests for the function
source("FunctionsLR.R")

# compatibility check
# all 1 col
X = matrix(c(rep(1, 5), rep(2, 5), rep(3, 5)), nrow = 5, ncol = 3)
any(X[ , 1] != 1)
X = matrix(c(rep(2, 5), rep(1, 5), rep(3, 5)), nrow = 5, ncol = 3)
any(X[ , 1] != 1)

# getting prob's from mat p
X = matrix(c(1,3,2,4,8,6,3,10,1,5,6,2), nrow = 4, ncol = 3)
X
y = c(0, 2, 0, 1)
# using techniques to access double array by addresses in C
t(X)[seq(0, nrow(X) - 1) * ncol(X) + (y + 1)]
seq(0, (nrow(X) - 1) * ncol(X), by = ncol(X))
seq(0, nrow(X) - 1) * ncol(X)
library(microbenchmark)
microbenchmark(
  seq(0, (nrow(X) - 1) * ncol(X), by = ncol(X)),
  seq(0, nrow(X) - 1) * ncol(X)
)
X[(seq(nrow(X)) - 1) * ncol(X) + y]

# creating indicator matrix
# indi_ij == indicator(y_i = j) where i:nsample and j:class
n = 10
k = 3
y = c(0,1,2,1,0,0,1,2,0,1)
y_adjusted = y + 1
indi = matrix(0, nrow = n, ncol = k)
indi[cbind(seq_along(y_adjusted), y_adjusted)] <- 1
length(y_adjusted)
indi

# Can we vectorize hessian? diag %*% mat
# diag = n * n
n = 4
p = 2
k = 3
# p mat is n * k
p_mat = matrix(seq(1,12), nrow = n, ncol = k)
p_mat
# d = seq(n)
# X = n * p
X = matrix(seq(1,8), nrow = n, ncol = p)
X
dim(outer(p_mat, X, "*"))[1,1, , ]
outer(p_mat, X, "*") [2, 1, , ]

mat_to_prod = matrix(1, nrow = 3, ncol = 2)
c(1,3,2) * mat_to_prod
crossprod(mat_to_prod, c(1,3,2) * mat_to_prod)
c(1,3,2) * (1 - c(2,1,3))

p_k_col = seq(n)
p_k_col
(p_k_col * (1-p_k_col)) * X
X

# test solve
solve(diag(seq(n), nrow = n), diag(1, nrow = n))

# tests on random matrix
n = 100
p = 10
k = 20
X_not_1 = matrix(rnorm(n * p), nrow = n, ncol = p)
y = sample(k, n, replace = TRUE) - 1
n_t = 20
X_t_not_1 = matrix(rnorm(n_t * p), nrow = n_t, ncol = p)
y_t = sample(k, n_t, replace = TRUE) - 1
# first row stop condition works?
LRMultiClass(X_not_1, y, X_t_not_1, y_t)
# Incompatible dimension
# generating function
X_generator <- function(n, p) {
  X = matrix(rnorm(n * (p - 1)), nrow = n, ncol = (p - 1))
  # 1's for first row
  X = cbind(rep(1, n), X)
  return(X)
}

# sample mismatch
missample_X = X_generator(n, p)
missample_Xt = X_generator(n_t + 1, p)
LRMultiClass(missample_X, y, missample_Xt, y_t)

# covariates mismatch
lackcovariate_X = X_generator(n, p - 1)
Xt = X_generator(n_t, p)
LRMultiClass(lackcovariate_X, y, Xt, y_t)

# normal cases
source("FunctionsLR.R")
X = X_generator(n, p)
Xt = X_generator(n_t, p)
# dim(X)
# dim(Xt)
LRMultiClass(X, y, Xt, y_t)

# debug
indi = matrix(0, nrow = n, ncol = k)
dim(indi)
print(length(y_adjusted))
print(seq_along(y_adjusted), y_adjusted)
indi[cbind(seq_along(y_adjusted), y_adjusted)] <- 1

# actual dataset: iris
iris_data <- iris
# convert y to 0~k class
iris_data$Species <- as.numeric(iris_data$Species) - 1
# permute
iris_data = iris_data[sample(nrow(iris_data)), ]
iris_data
# divide train and test
iris_train = iris_data[1:((0.8) * nrow(iris_data)), ]
iris_train
iris_test = iris_data[(0.8 * nrow(iris_data)):nrow(iris_data), ]
iris_test
# y
iris_y = as.vector(iris_train$Species)
iris_yt = as.vector(iris_test$Species)
# X
iris_X = cbind(rep(1, nrow(iris_train)), iris_train[, -ncol(iris_train)])
iris_X = as.matrix(iris_X)
iris_Xt = cbind(rep(1, nrow(iris_test)), iris_test[, -ncol(iris_test)])
iris_Xt = as.matrix(iris_Xt)

# compare performance
library(nnet)
iris_package = cbind(rep(1, nrow(iris_train)), iris_train)
# iris_package
package_model <- multinom(Species ~ ., data = iris_package)
my_model <- LRMultiClass(iris_X, iris_y, iris_Xt, iris_yt, numIter = 100)
# error by package
(sum(max.col(package_model$fitted.values) != (iris_y + 1)) / length(iris_y)) * 100
# error of mine
my_model
my_model$beta
summary(package_model)
microbenchmark(
  LRMultiClass(iris_X, iris_y, iris_Xt, iris_yt, numIter = 100), 
  multinom(Species ~ ., data = iris_package)
)