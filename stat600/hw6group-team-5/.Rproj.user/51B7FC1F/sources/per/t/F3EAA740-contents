# Application of K-means algorithm to ZIPCODE data
source("FunctionsKmeans.R")
# Rand Index Calculation example
require(fossil)
cluster1 <- c(2,2,1,3)
cluster2 <- c(3,3,2,1)
rand.index(cluster1, cluster2) # clusters match

# Load the ZIPCODE data
zipcode <- read.table("ZIPCODE.txt", header = F)

# Extract the true digits
Y1 <- zipcode[ , 1]

# Extract the data points
X1 <- zipcode[ , -1]
n1 = nrow(X1)
p1 = ncol(X1)
K1 = 10

set.seed(1)
res_base = MyKmeans(X1, K1)
res_cur = GroupHW::MyKmeans(X1, K1)
table(res_base)
table(res_cur)
# [ToDo] Try K-means algorithm nRep times with different starting points on ZIPCODE data. Calculate Rand Index at each replication
nRep <- 5
rand_index = 0
cnt = 0
# iterate for 50 times
for(i in 1:nRep) {
  trial <- try(y_hat <<- MyKmeans(X1, K1))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt = cnt + 1
  rand_index = rand_index + rand.index(Y1, y_hat)
}
rand_index = rand_index / cnt
rand_index
cnt

# library(GroupHW)
rand_index_c = 0
cnt_c = 0
for(i in 1:nRep) {
  trial <- try(y_hat_c <<- GroupHW::MyKmeans(X1, K1))
  # skip the failed one
  if(inherits(trial, "try-error")) {
    next
  }
  cnt_c = cnt_c + 1
  rand_index_c = rand_index_c + rand.index(Y1, y_hat_c)
}

# [ToDo] Report mean Rand Index
rand_index_c = rand_index_c / cnt_c
rand_index_c
cnt_c

## results
# > rand_index
# [1] 0.9105936
# > cnt
# [1] 50

### Also, I also ran comparisons with base kmeans, and it's also in tests.r

# [ToDo] Report mean run time for one replication on your machine
library(microbenchmark)
microbenchmark(
  MyKmeans(X1, K1),
  GroupHW::MyKmeans(X1, K1),
  times = 5
)
#         expr                                                  min       lq     mean
# MyKmeans(X, K)                                              1.051792 1.750646 1.913303
# kmeans(X, centers = 10, iter.max = 100, algorithm = "Lloyd") 1.083030 1.372673 1.666200
# median       uq      max neval
# 1.891406 2.120334 2.780963    10
# 1.483034 1.864435 3.156567    10