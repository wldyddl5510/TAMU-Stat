# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")
set.seed(1)
# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
rib_fitlasso = fitLASSO(X, Y, NULL, n_lambda = 60)

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
num_of_nonzero = colSums(rib_fitlasso$beta_mat != 0)
plot(rib_fitlasso$lambda_seq, num_of_nonzero)

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark)
microbenchmark(
  fitLASSO(X, Y, NULL, n_lambda = 60),
  times = 10
)

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
# Unit: seconds
# expr      min       lq     mean   median       uq      max
# 2.959695 3.114957 3.634137 3.393925 3.495275 6.587223
# neval
# 10

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
rib_cvlasso = cvLASSO(X, Y, NULL, 30)
rib_cvlasso

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
plot(rib_cvlasso$lambda_seq, rib_cvlasso$cvm, col = "red", ylim = c(0,1))
lines(rib_cvlasso$lambda_seq, rib_cvlasso$cvm + rib_cvlasso$cvse)
lines(rib_cvlasso$lambda_seq, rib_cvlasso$cvm - rib_cvlasso$cvse)
rib_cvlasso$lambda_min
rib_cvlasso$lambda_1se
