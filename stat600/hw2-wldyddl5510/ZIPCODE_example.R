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
Y <- zipcode[ , 1]

# Extract the data points
X <- zipcode[ , -1]
n = nrow(X)
p = ncol(X)
K = 10

# [ToDo] Try K-means algorithm nRep times with different starting points on ZIPCODE data. Calculate Rand Index at each replication
nRep <- 50
rand_index = 0
cnt = 0
# iterate for 50 times
for(i in 1:nRep) {
  trial <- try(y_hat <<- MyKmeans(X, K))
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

## results
# > rand_index
# [1] 0.9105936
# > cnt
# [1] 50

### Also, I also ran comparisons with base kmeans, and it's also in tests.r

# [ToDo] Report mean run time for one replication on your machine
library(microbenchmark)
microbenchmark(
  MyKmeans(X, K),
  # kmeans(X, centers = 10, iter.max = 100, algorithm = "Lloyd"),
  times = 10
)
#         expr                                                  min       lq     mean
# MyKmeans(X, K)                                              1.051792 1.750646 1.913303
# kmeans(X, centers = 10, iter.max = 100, algorithm = "Lloyd") 1.083030 1.372673 1.666200
# median       uq      max neval
# 1.891406 2.120334 2.780963    10
# 1.483034 1.864435 3.156567    10