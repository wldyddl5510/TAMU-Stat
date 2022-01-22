# Code to apply robust PCA via ADMM to a given input

###################################
## Generate the simulated data
###################################
# Create matrix M according to M = low rank + few large elements
n = 100
p = 40

# Create rank 3 component
set.seed(1234)
out <- svd(matrix(rnorm(n*p), n, p), nu = 3, nv = 3)
trueL <- out$u %*% diag(out$d[1:3]) %*% t(out$v)

# Create sparse component
trueS <- matrix(rt(n*p, df = 1), n, p)
trueS[abs(trueS) < 2.5] <- 0
sum(trueS !=0)

# Create M by combining the above
M <- trueL + trueS
# Verify that M itself is not small rank and not even close to being rank 3
svd(M)$d

###################################
## Apply ADMM
###################################
source("ADMMfunctions.R")

outADMM <- robustPCAadmm(M, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1, eps = 0.0001)

# What is the rank of the solution?
round(svd(outADMM$L)$d, 2)
# How many non-zero components are in the sparse part?
sum(outADMM$S != 0)
# Is the constraint satisfied approximately?
sum((M - outADMM$S - outADMM$L)^2)

# Check the code speed (only after verified the code is correct)
library(microbenchmark)
microbenchmark(
  robustPCAadmm(M, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1, eps = 0.0001),
  times = 10
)

# Unit: seconds
# expr
# robustPCAadmm(M, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1,      eps = 1e-04)
# min       lq     mean   median       uq      max neval
# 37.7989 38.46728 40.81174 38.91679 43.38914 47.69203    10