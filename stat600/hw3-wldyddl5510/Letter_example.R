# Application of multi-class logistic to letters data

# Load the letter data
#########################
# Training data
letter_train <- read.table("Data/letter-train.txt", header = F, colClasses = "numeric")
Y <- letter_train[, 1]
X <- as.matrix(letter_train[, -1])

# Testing data
letter_test <- read.table("Data/letter-test.txt", header = F, colClasses = "numeric")
Yt <- letter_test[, 1]
Xt <- as.matrix(letter_test[, -1])

# [ToDo] Make sure to add column for an intercept to X and Xt
n = nrow(X)
ntest = nrow(Xt)
X = cbind(rep(1, n), X)
Xt = cbind(rep(1, ntest), Xt)

# Source the LR function
source("FunctionsLR.R")

# [ToDo] Try the algorithm LRMultiClass with lambda = 1 and 50 iterations. Call the resulting object out, i.e. out <- LRMultiClass(...)
out = LRMultiClass(X, Y, Xt, Yt)

# The code below will draw pictures of objective function, as well as train/test error over the iterations
plot(out$objective, type = 'o')
plot(out$error_train, type = 'o')
plot(out$error_test, type = 'o')
# result of the last step
print(out$objective[length(out$objective)])
print(out$error_train[length(out$error_train)])
print(out$error_test[length(out$error_test)])

# Feel free to modify the code above for different lambda/eta/numIter values to see how it affects the convergence as well as train/test errors
# larger and smaller lambda
out_lambda2 = LRMultiClass(X, Y, Xt, Yt, lambda = 2)
out_lambda05 = LRMultiClass(X, Y, Xt, Yt, lambda = 0.5)

# objective with diff lambda
plot(out_lambda2$objective, type = 'o')
plot(out_lambda05$objective, type = 'o')
# compare last step
print(out_lambda2$objective[length(out_lambda2$objective)])
print(out_lambda05$objective[length(out_lambda05$objective)])

# train error with diff lambda
plot(out_lambda2$error_train, type = 'o')
plot(out_lambda05$error_train, type = 'o')
print(out_lambda2$error_train[length(out_lambda2$error_train)])
print(out_lambda05$error_train[length(out_lambda05$error_train)])

# test error with diff lambda
plot(out_lambda2$error_test, type = 'o')
plot(out_lambda05$error_test, type = 'o')
print(out_lambda2$error_test[length(out_lambda2$error_test)])
print(out_lambda05$error_test[length(out_lambda05$error_test)])

# larger and smaller eta
out_eta05 = LRMultiClass(X, Y, Xt, Yt, eta = 0.5)
out_eta001 = LRMultiClass(X, Y, Xt, Yt, eta = 0.01)

# objective with diff eta
plot(out_eta05$objective, type = 'o')
plot(out_eta001$objective, type = 'o')
print(out_eta05$objective[length(out_eta05$objective)])
print(out_eta001$objective[length(out_eta001$objective)])

# train error with diff eta
plot(out_eta05$error_train, type = 'o')
plot(out_eta001$error_train, type = 'o')
print(out_eta05$error_train[length(out_eta05$error_train)])
print(out_eta001$error_train[length(out_eta001$error_train)])

# test error with diff eta
plot(out_eta05$error_test, type = 'o')
plot(out_eta001$error_test, type = 'o')
print(out_eta05$error_test[length(out_eta05$error_test)])
print(out_eta001$error_test[length(out_eta001$error_test)])

# larger and smaller Numiter
out_iter25 = LRMultiClass(X, Y, Xt, Yt, numIter = 25)
out_iter100 = LRMultiClass(X, Y, Xt, Yt, numIter = 100)

# objective with diff numIter
plot(out_iter25$objective, type = 'o')
plot(out_iter100$objective, type = 'o')
print(out_iter25$objective[length(out_iter25$objective)])
print(out_iter100$objective[length(out_iter100$objective)])

# train error with diff iter
plot(out_iter25$error_train, type = 'o')
plot(out_iter100$error_train, type = 'o')
print(out_iter25$error_train[length(out_iter25$error_train)])
print(out_iter100$error_train[length(out_iter100$error_train)])

# test error with diff iter
plot(out_iter25$error_test, type = 'o')
plot(out_iter100$error_test, type = 'o')
print(out_iter25$error_test[length(out_iter25$error_test)])
print(out_iter100$error_test[length(out_iter100$error_test)])

# [ToDo] Use microbenchmark to time your code with lambda=1 and 50 iterations. To save time, only apply microbenchmark 5 times.
library(microbenchmark)
microbenchmark(
  LRMultiClass(X, Y, Xt, Yt, lambda = 1, numIter = 50),
  svd(matrix(rnorm(50000), 1000, 500)), 
  times = 5
)

# [ToDo] Report the median time of your code from microbenchmark above in the comments below

# Median time:  3342.1790 (in milliseconds)