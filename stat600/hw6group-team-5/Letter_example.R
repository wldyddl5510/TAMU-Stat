# Application of multi-class logistic to letters data

# Load the letter data
#########################
# Training data
letter_train <- read.table("letter-train.txt", header = F, colClasses = "numeric")
Y1 <- letter_train[, 1]
X1 <- as.matrix(letter_train[, -1])

# [ToDo] Make sure to add column for an intercept to X and Xt
n1 = nrow(X1)
X1 = cbind(rep(1, n1), X1)

# Source the LR function
source("FunctionsLR.R")

# [ToDo] Try the algorithm LRMultiClass with lambda = 1 and 50 iterations. Call the resulting object out, i.e. out <- LRMultiClass(...)
out = LRMultiClass(X1, Y1)

out_c = GroupHW::LRMultiClass(X1, Y1)
# The code below will draw pictures of objective function, as well as train/test error over the iterations
plot(out$objective, type = 'o')
plot(out_c$objective, type = 'o')
# result of the last step
all.equal(out$objective, as.numeric(out_c$objective))
all.equal(out$beta, out_c$beta)
print(out$objective[length(out$objective)])
print(out_c$objective[length(out$objective)])

# Feel free to modify the code above for different lambda/eta/numIter values to see how it affects the convergence as well as train/test errors
# larger and smaller lambda
out_lambda2 = LRMultiClass(X1, Y1, lambda = 2)
out_lambda05 = LRMultiClass(X1, Y1, lambda = 0.5)
out_c_lambda2 = GroupHW::LRMultiClass(X1, Y1, lambda = 2)
out_c_lambda05 = GroupHW::LRMultiClass(X1, Y1, lambda = 0.5)

# objective with diff lambda
plot(out_lambda2$objective, type = 'o')
plot(out_c_lambda2$objective, type = 'o')
plot(out_lambda05$objective, type = 'o')
plot(out_c_lambda05$objective, type = 'o')
# compare last step
all.equal(out_lambda2$objective, as.numeric(out_c_lambda2$objective))
all.equal(out_lambda2$beta, out_c_lambda2$beta)
all.equal(out_lambda05$objective, as.numeric(out_c_lambda05$objective))
all.equal(out_lambda05$beta, out_c_lambda05$beta)

# train error with diff lambda
#plot(out_lambda2$error_train, type = 'o')
#plot(out_lambda05$error_train, type = 'o')
#print(out_lambda2$error_train[length(out_lambda2$error_train)])
#print(out_lambda05$error_train[length(out_lambda05$error_train)])

# test error with diff lambda
#plot(out_lambda2$error_test, type = 'o')
#plot(out_lambda05$error_test, type = 'o')
#print(out_lambda2$error_test[length(out_lambda2$error_test)])
#print(out_lambda05$error_test[length(out_lambda05$error_test)])

# larger and smaller eta
out_eta05 = LRMultiClass(X1, Y1, eta = 0.5)
out_c_eta05 = GroupHW::LRMultiClass(X1, Y1, eta = 0.5)
out_eta001 = LRMultiClass(X1, Y1, eta = 0.01)
out_c_eta001 = GroupHW::LRMultiClass(X1, Y1, eta = 0.01)
all.equal(out_eta05$objective, as.numeric(out_c_eta05$objective))
all.equal(out_eta001$objective, as.numeric(out_c_eta001$objective))
all.equal(out_eta05$beta, out_c_eta05$beta)
all.equal(out_eta001$beta, out_c_eta001$beta)


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
out_iter25 = LRMultiClass(X1, Y1, numIter = 25)
out_iter100 = LRMultiClass(X1, Y1, numIter = 100)
out_iter25_c = GroupHW::LRMultiClass(X1, Y1, numIter = 25)
out_iter100_c = GroupHW::LRMultiClass(X1, Y1, numIter = 100)
all.equal(out_iter25$beta, out_iter25_c$beta)
all.equal(out_iter25$objective, as.numeric(out_iter25_c$objective))
all.equal(out_iter100$objective, as.numeric(out_iter100_c$objective))
all.equal(out_iter100$beta, out_iter100_c$beta)

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
  LRMultiClass(X1, Y1, lambda = 1, numIter = 50),
  GroupHW::LRMultiClass(X1, Y1, lambda = 1, numIter = 50),
  times = 10
)

# [ToDo] Report the median time of your code from microbenchmark above in the comments below

# Median time:  3342.1790 (in milliseconds)