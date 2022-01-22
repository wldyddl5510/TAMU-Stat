# Generate data from linear regression model and calculate the least squares vector of coefficients
#####################################################################################################

# Source your functions
source("FunctionsLM.R")

# Model parameters
p = 10 # number of covariates
sigma = 2 # noise standard deviation
beta = c(2, rep(0, p-1)) # true vector of coefficients

# Training data generator
n = 100 # sample size for training data
X = matrix(rnorm(n * p), n, p) # n by p matrix of predictors
# [ToDo] Use generateY function to generate Y for training data with default seed
Y = generateY(X, beta, sigma)

# [ToDo] Use calculateBeta function to calculate beta_LS
beta_LS = calculateBeta(X, Y)

# [ToDo] Use calculateEstimationError to assess the estimation error measured by squared eucledian distance - ||beta - beta_LS||_2^2. Report the error in the comments.
beta_L2_error = calculateEstimationError(beta, beta_LS)
print(beta_L2_error)
# beta_L2_error == 1.009905

# Testing data generator
n = 200 # sample size for testing data
Xtest = matrix(rnorm(n * p), n, p) # n by p matrix of covariates
# [ToDo] Use generateY function to generate Ytest for testing data with seed = 678910
Ytest = generateY(Xtest, beta, sigma, seed = 678910)

# [ToDo] Use calculatePredictionError to asses the prediction error on Ytest. Report the error in the comments.
# calculate the beta_LS
beta_LS_test = calculateBeta(Xtest, Ytest)
# calculate the pred error for given X and Beta
pred_error_test = calculatePredictionError(Ytest, Xtest, beta_LS_test)
print(pred_error_test)
# pred_error_test == 713.432

# [ToDo] Use calculatePredictionError to asses the prediction error on Ytest based only on the first covariate. Report the error in the comments.
# Hint: to avoid error of non-conformable arguments, use Xtest[, 1, drop = FALSE]
# we use Xtest[, 1, drop = FALSE] to extract 1st cov
Xtest_first_cov = Xtest[, 1, drop = FALSE]
# calculate the beta for only first cov
beta_LS_test_first_cov = calculateBeta(Xtest_first_cov, Ytest)
# calculate the pred error for only first cov
pred_error_test_first_cov = calculatePredictionError(Ytest, Xtest_first_cov, beta_LS_test_first_cov)
print(pred_error_test_first_cov)
# pred_error_test_first_cov == 744.4918