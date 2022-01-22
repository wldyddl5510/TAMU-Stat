# Matrix soft-thresholding (applying soft-thresholding to each element of a matrix)
# INPUT
# A - n x p matrix
# lambda - non-negative thresholding parameter
# OUTPUT
# n x p matrix, each element is soft-thresholded element of A matrix
soft <- function(A, lambda){
  # Fill in
  # sgn(x) * max(|x| - lambda, 0) all element-wise
  return(sign(A) * pmax(abs(A) - lambda, 0))
}

# Proximal operator for nuclear norm (matrix with soft-thresholded singular values)
# INPUT
# A - n x p matrix
# lambda - non-negative thresholding parameter
# OUTPUT
# A list with
#  newA - n x p matrix, same singular vectors as A, but singular values soft-thresholded using lambda
#  newd - a vector of soft-thresholded singular values, that is singular values of newA (returning it here saves one SVD step later)
soft_nuclear <- function(A, lambda){
  # Fill in
  # compute svd
  svd_of_A = svd(A)
  # soft(D) = newD
  newd = soft(svd_of_A$d, lambda)
  # U * newD * V^t = newA
  u = svd_of_A$u
  newA = tcrossprod(u * rep(newd, rep(nrow(u), ncol(u))), svd_of_A$v)
  # Returns a list
  # newd - soft-thresholded singular values, vector
  # newA - n x p matrix with singular values newd
  return(list(newd = newd, newA = newA))
}

# Objective function value of robustPCA, |L|_* + gamma|S|_1,
# INPUT
# Ld - a vector of singular values of L (this avoids recalculating svd more than needed)
# S - a n by p matrix, current value of S
# gamma - a positive scalar, current value of gamma
# OUTPUT
# A scalar corresponding to |L|_* + gamma|S|_1
robustPCAobj <- function(Ld, S, gamma = 0.1){
  # Fill in
  # nuclear norm
  L_nuclear_norm = sum(Ld)
  # L1 norm
  S_one_norm = sum(abs(S))
  # sum
  return(L_nuclear_norm + gamma * S_one_norm)
}


# Scaled ADMM algorithm for Robust PCA
# INPUT
# M - a n x p matrix
# gamma - positive scalar, default value 0.1
# Sinit - a n x p matrix, starting value for S. If none is supplied, initialize with matrix of zeros. If supplied, check for compatibility of dimensions with M.
# etainit - a n x p matrix, starting value of eta. If none is supplied, initialize with matrix of zeros. If supplied, check for compatibility of dimensions with M.
# tau - a positive scalar, ADMM parameter for scaled ADMM version, default value is 1.
# eps - a positive scalar, convergence tolerance criteria (difference in objective function values), default value is 0.001
# OUTPUT
# A list with
#   L - value of matrix L at convergence
#   S - value of matrix S at congergence
#   eta - value of matrix eta at convergence
robustPCAadmm <- function(M, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1, eps = 0.001){
  # Check for compatibility between supplied Sinit, etainit and M. If Sinit or etainit is NULL, initialize the corresponding matrix as a matrix of zeros.
  n = nrow(M)
  p = ncol(M)
  # null
  if(is.null(Sinit)) {
    Sinit = matrix(0, n, p)
  } else { # compatibility check
    if(nrow(Sinit) != n) {
      stop("Num of rows between M and S must match")
    }
    if(ncol(Sinit) != p){
      stop("Num of columns between M and S must match")
    }
  }
  if(is.null(etainit)) {
    etainit = matrix(0, n, p)
  } else { # compatibility check
    if(nrow(etainit) != n) {
      stop("Num of rows between M and eta must match")
    }
    if(ncol(etainit) != p){
      stop("Num of columns between M and eta must match")
    }
  }
  # Initialize L as L = M - S
  S = Sinit
  eta = etainit
  L = M - S
  # calculate Ld for current L
  # calculate init objective
  curr_obj = robustPCAobj(svd(L)$d, S, gamma)
  diff = curr_obj
  ## Implement ADMM algorithm, which alternates updates of L, S and eta until convergence (first change in objective function values less than eps)
  while(diff >= eps) {
    # update L
    L_d_A = soft_nuclear(M - S - eta, tau)
    L = L_d_A$newA
    Ld = L_d_A$newd
    # update S
    S = soft(M - L - eta, gamma * tau)
    # update eta
    eta = eta + (L + S - M) / tau
    # update objective
    new_obj = robustPCAobj(Ld, S, gamma)
    # calculate the difference
    diff = abs(new_obj - curr_obj)
    curr_obj = new_obj
  }
  # Return L, S and eta at convergence
  return(list(L = L, S = S, eta = eta))
}