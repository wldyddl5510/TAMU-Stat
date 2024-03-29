# Use this file to create tests/debug your functions
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
# Source the functions
source("FunctionsKmeans.R")
source("kmeans_c_wrapper.R")
library(fossil)

# packages
library(microbenchmark)

# comparison for distance matrix computation: X1 - M
n = 4
p = 3
k = 2
X1 = matrix(c(1:(n * p)), n, p)
M1 = matrix(c(1:(k*p)), k, p)
# M1 = matrix(c((k * p): 1), p, k)
# M1 / rep(c(1,3), each = nrow(M1))
diff_mat = M1 - rep(X1[1, ], each = nrow(M1))
diff_mat^2
# method 1
which.min(diag(tcrossprod(M1 - rep(X1[1, ], each = nrow(M1)))))
# method 2
rowSums(diff_mat^2)
# compare
microbenchmark(
  which.min(diag(tcrossprod(diff_mat))),
  which.min(rowSums(diff_mat^2))
)

# taking min
microbenchmark(
  which.min(rowSums(diff_mat^2)),
  max.col(-rowSums(diff_mat^2))
)

# ver2 of calculating dist
L2_dist <- function(v, w){
  return(sum((v - w)^2))
}

# each element in (i,j) will stand for dist(X_i, M_j)
mat_dist1 <- function(X, M){
  dist_mat <- apply(M, 1, function(y) apply(X, 1, L2_dist, y))
  return(dist_mat)
}

mat_dist2 <- function(X, M) {
  dist_mat <- outer(rowSums(X^2), rowSums(M^2), '+') - tcrossprod(X, 2 * M)
  return(dist_mat)
}

mat_dist2(X1, M1)
rowSums(X1^2)
rowSums(M1^2)
tcrossprod(X1, 2 * M1)
M1
# calculate maxcol of -dist_mat
whole_min_cluster1 <- function(X, M) {
  return(max.col(-mat_dist1(X, M)))
}

whole_min_cluster2 <- function(X, M) {
  return(max.col(-mat_dist2(X, M)))
}

whole_min_cluster2(X1, M1)

# calculating
mat_dist1(X1, M1)
X1
M1
sum((X1[2,]-M1[2,])^2)
mat_dist2(X1, M1)

max.col(-mat_dist2(X1, M1))

#indexing
A = mat_dist2(X1, M1)
max.col(-A)
colMeans(X1[which(max.col(-A) %in% 1), ])


# rowwise iteration version
# where is each row
# for(X_row_i in 1:X_nrow) {
#   # select the row
#    i_row = X[X_row_i, ]
#    # get the cluster of given row
#    # calculate the diff
#    # diff_mat = (M - rep(i_row, each = K))
#    # i_row_cluster = which.min(rowSums(diff_mat^2))
#    i_row_cluster = row_min_cluster(i_row, M, K)
#    # fil in the result vector
#    Y[X_row_i] = i_row_cluster
#    # kth cluster exists
#    cluster_indicator[i_row_cluster] = TRUE
#    # summing its value to its cluster
#    cluster_sum[i_row_cluster, ] <- cluster_sum[i_row_cluster, ] + i_row
#    cluster_cnt[i_row_cluster] <- cluster_cnt[i_row_cluster] + 1
# }

# cluster disappears
# # if(any(cluster_indicator) == FALSE) {
#   stop("certain cluster has been disappeared")
# }
# # reassign cluster mean
# new_mean <- (cluster_sum / cluster_cnt)
# 

# double for loop version
#   X_row = X[X_row_i, ]
#   dist_between_each_clus = rep(NA, K)
#   for (M_row_j in 1:K) {
#     M_row = M[M_row_j, ]
#     # print(dist_between_each_clus)
#     dist_between_each_clus[M_row_j] = dist_row_by_loop(X_row, M_row)
#   }
#   # X_row_cluster = max.col(-dist_between_each_clus)
#   X_row_cluster = which.min(dist_between_each_clus)
#   cluster_indicator[X_row_cluster] = TRUE
#   Y[X_row_i] = X_row_cluster
#   cluster_sum[X_row_cluster, ] <- cluster_sum[X_row_cluster, ] + X_row
#   cluster_cnt[X_row_cluster] <- cluster_cnt[X_row_cluster] + 1
# }

# all vs any
m = 1000000
boolean_true = rep(TRUE, m)
boolean_false = rep(FALSE, m)
microbenchmark(
  all(boolean_true),
  !any(boolean_false)
)
nrow(NULL)
ncol(NULL)

!all(c(, TRUE))

# all equal vs identical
p = 400
k = 30
X = matrix(rep(0, p * k), k, p)
Y = matrix(rep(0, p * k), k, p)
# which[X == Y]

# actual example
n = 40
p = 20
K = 10
X = matrix(rnorm(n * p), n, p)
M = matrix(rnorm(K * p), K, p)

all.equal(MyKmeans(X, K, M), as.numeric(MyKmeans_c(X, K, M)))
sourceCpp("kmeanscpp.cpp")
all.equal(as.numeric(test_rows(X, c(0, 1))), colMeans(X[1:2, ]))
X[2, ]
set.seed(1)
whole_min_cluster(X, M)
as.numeric(whole_min_cluster_c(X, M, n, K))
all.equal(whole_min_cluster(X, M),as.numeric(whole_min_cluster_c(X, M, n, K)))
M_invalid_row = matrix(rnorm((K - 1) * p), K - 1, p)
M_invalid_col = matrix(rnorm(K * (p - 1)), K, p - 1)
M_NULL = NULL
# compare with Kmeans in base
# error cases
MyKmeans(X, K, M_invalid_row)
MyKmeans(X, K, M_invalid_col)
# NULL cases
MyKmeans(X, K, M_NULL)
# small case
MyKmeans(X1, 2, M1)
# K == 1
MyKmeans(X, 1, M)
Rprof(gc.profiling = TRUE)
MyKmeans(X, K, M)
MyKmeans_wrapper(X, K, M)
Rprof(NULL) # stop monitoring
summaryRprof() # see the report
# mean(MyKmeans(X, K, M) == kmeans(X, M, iter.max = 100, algorithm = "Lloyd")$cluster)
rand.index(MyKmeans(X, K, M), kmeans(X, M, iter.max = 100, algorithm = "Lloyd")$cluster)
library(microbenchmark)
microbenchmark(
  MyKmeans(X, K, M),
  kmeans(X, M, iter.max = 100, algorithm = "Lloyd")$cluster,
  MyKmeans_wrapper(X, K, M),
  times = 5
)

#####ZIPCODE tests
# Load the ZIPCODE data
zipcode <- read.table("ZIPCODE.txt", header = F)

# Extract the true digits
Y <- zipcode[ , 1]

# Extract the data points
X <- zipcode[ , -1]
n = nrow(X)
p = ncol(X)
K = 10

nRep <- 10
rand_index = 0
cnt = 0

# iterate for 10 times
for(i in 1:nRep) {
  trial <- try(y_hat <<- MyKmeans(X, K))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt = cnt + 1
  rand_index = rand_index + rand.index(Y, y_hat)
}

rand_index = rand_index / cnt
rand_index
cnt

rand_index_c = 0
cnt_c = 0
# iterate for 10 times
for(i in 1:nRep) {
  trial <- try(y_hat <<- MyKmeans_wrapper(X, K))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt_c = cnt_c + 1
  rand_index_c = rand_index_c + rand.index(Y, y_hat)
}

rand_index_c = rand_index_c / cnt_c
rand_index_c
cnt_c

# iterate for 50 times
for(i in 1:nRep) {
  M = matrix(rnorm(K * p), K, p)
  trial <- try(y_hat <<- MyKmeans(X, K, M))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt = cnt + 1
  rand_index = rand_index + rand.index(Y, y_hat)
}

# [ToDo] Report mean Rand Index
rand_index = rand_index / cnt
rand_index
cnt

### I tried M = NULL with 50 iterations here.
nRep <- 50
rand_index_null = 0
cnt_null = 0
# iterate for 50 times
for(i in 1:nRep) {
  trial <- try(y_hat <<- MyKmeans(X, K))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt_null = cnt_null + 1
  rand_index_null = rand_index_null + rand.index(Y, y_hat)
}

rand_index_null = (rand_index_null / cnt_null)
rand_index_null
cnt_null

# > rand_index_null
# [1] 0.9099871
# > cnt_null
# [1] 50

### Also, I also ran comparison with base functions
cnt_base = 0
rand_index_base = 0
for(i in 1:nRep) {
  M = matrix(rnorm(K * p), K, p)
  y_hat_base = kmeans(X, M, iter.max = 100, algorithm = "Lloyd")$cluster
  cnt_base = cnt_base + 1
  rand_index_base = rand_index_base + rand.index(Y, y_hat_base)
}
rand_index_base = (rand_index_base / cnt_base)
rand_index_base
cnt_base

# > rand_index_base
# [1] 0.9049956
# > cnt_base
# [1] 50

## base with init on NULL
cnt_base = 0
rand_index_base = 0
for(i in 1:nRep) {
  y_hat_base = kmeans(X, centers = 10, iter.max = 100, algorithm = "Lloyd")$cluster
  cnt_base = cnt_base + 1
  rand_index_base = rand_index_base + rand.index(Y, y_hat_base)
}
rand_index_base = (rand_index_base / cnt_base)
rand_index_base
cnt_base

# > rand_index_base
# [1] 0.9111734
# > cnt_base
# [1] 50