## Homework 5 - LASSO implementation in Rcpp with Armadillo library

# Introduction
In this homework, you will practice Rcpp and the use of Armadillo library by implementing the LASSO coordinate-descent algorithm from HW4. We will only focus on `fitLASSOstandardized` functions and corresponding helpers as you will see the biggest speed improvement there. For simplicity, we will also **avoid doing any compatibility/input checks** in C++. In practice, you will often have an R code wrapper with all compatibility checks that then calls the corresponding function in C++.

## Starter code

**LassoInC.cpp** contains the starter for C++ code. You will have to modify this starter to create the following functions:

  - `soft_c` - this is the C++ version of the `soft` function from HW4. It should return the soft-thresholded value of a. 
  
  - `lasso_c` - this is the C++ version of the `lasso` function from HW4. It should return the value of lasso objective function at given arguments.
  
  - `fitLASSOstandardized_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from HW4. The  modification: it only returns $\beta$ vector (it no longer returns objective function value)
  
  - `fitLASSOstandardized_seq_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from HW4. The modification: it only takes used-supplied lambda sequence and it only returns matrix of $\beta$ values. You can assume that the user-supplied sequence has already been sorted from largest to smallest.
  
All these functions should be directly available from R after executing **sourceCpp** command with corresponding source file (see **TestsCpp.R**). You are welcome to create additional C++ functions within the same file if you want, we will only explicitly test the above 4 functions. You can assume that all the inputs to these functions are correct (e.g. dimensions match, no negative tuning parameters, etc.)

**TestsCpp.R**  contains starter file for testing your C++ functions against your R functions. You should do the following

  - upload your **LASSOfunctions.R** script from HW4 to this project repository
  - develop at least 2 tests for each function **to check equality of the corresponding outputs between your R and C++ versions**
  - do 2 microbenchmarks comparisons as per comments
  - do speed test on riboflavin data as coded in the end of that file
  
Things to keep in mind when implementing:

 - you will use Armadillo library classes for matrices and vectors, and because of this, the whole Armadillo library of functions is at your disposal. As it's impossible to cover all functions in class, a part of this assignment is learning how to search for needed functions on your own using [Armadillo documentation](http://arma.sourceforge.net/docs.html). Using Ctrl(Command) + F on that page with a key word for the function you need, e.g. "inverse", will help you narrow down the function list for what you need.
 
 - one of the best ways to quickly learn how to translate the code is to study some relevant examples. In addition to examples in class, you may find the following [gallery of Rcpp and RcppArmadillo examples](https://gallery.rcpp.org) useful

 - test small chunks before testing larger chunks as it's harder to debug C++ code, and it's easier to miss a mistake in C++. I recommend first testing that a particular line of code does exactly what you think it does before going to the whole function
 
 - there are many more ways in C++ code to achieve a similar outcome than in R so we will be expecting to see much more differences across the assignments. If your code is heavily based on particular examples you found online, you will be expected to indicate this in the comments


  
## Grading for this assignment

Your assignment will be graded based on 

 * correctness _(50% of the grade)_

Take advantage of the fact that you have a (hopefully) correct R function to extensively test your C++ code. You also have to develop intput/output tests yourself for this particular assignment.
 
 * speed _(30% of the grade)_ 
 
You should expect to see a dramatic speed improvement when moving your code to C++. On riboflavin data, my C++ `fitLASSOstandardized_seq_c` code with 30 lambdas is around 40! times faster than the corresponding R version. You will get full points if your C++ code is **at most twice slower** comparable to my C++. You will loose 5 points for every fold over. You will get +5 bonus points if your **completely correct** code on 1st submission is faster than mine (median your time/median mine time < 0.9). 

  * code style/documentation _(10% of the grade)_

You need to comment different parts of the code so it's clear what they do, have good indentation, readable code with names that make sense. See guidelines on R style (similar guidelines apply to C++ code), and posted grading rubric.

* version control/commit practices _(10% of the grade)_
 
 I expect you to start early on this assignment, and work gradually. You want to commit often, have logically organized commits with short description that makes sense. See guidelines on good commit practices, and posted grading rubric.