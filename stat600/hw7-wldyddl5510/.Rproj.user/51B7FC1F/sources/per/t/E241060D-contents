# This is an empty script to store your own tests
source("ADMMfunctions.R")


#####################3
n1 = 50
p1 = 20
A = matrix(rnorm(n1 * p1), n1, p1)

res = svd(A)

library(microbenchmark)
microbenchmark(
  res$u %*% diag(res$d) %*% t(res$v),
  tcrossprod(res$u * rep(res$d, rep(nrow(res$u), ncol(res$u))), res$v),
  times = 5
)

all.equal(soft(diag(res$d), 1), diag(soft(res$d, 1)))
res$u %*% diag(res$d)
res$u * rep(res$d, rep(nrow(res$u), ncol(res$u)))
tcrossprod(res$u %*% diag(res$d), res$v)

all.equal(svd(soft_nuclear(A, 1)$newA)$d, soft_nuclear(A, 1)$newd)



#########################









n1 = 100
p1 = 10

set.seed(1)
n_rank = 7
out <- svd(matrix(rnorm(n1*p1), n1, p1), nu = 7, nv = 7)
out$d
trueL <- out$u %*% diag(out$d[1:7]) %*% t(out$v)

trueS <- matrix(rt(n1*p1, df = 1), n1, p1)
trueS[abs(trueS) < 2.5] <- 0
sum(trueS > 1e-50)

# Create M by combining the above
M1 <- trueL + trueS
# Verify that M itself is not small rank and not even close to being rank 7
svd(M1)$d

outADMM <- robustPCAadmm(M1, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1, eps = 0.0001)

# What is the rank of the solution?
round(svd(outADMM$L)$d, 2)
# How many non-zero components are in the sparse part?
sum(outADMM$S > 1e-50)
# Is the constraint satisfied approximately?
sum((M1 - outADMM$S - outADMM$L)^2)

