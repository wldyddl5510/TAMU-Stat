
# Homework 7 - ADMM algorithm to solve Robust PCA problem

**Attention:** Because math rendering of .Rmd is not ideal, please see
the enclosed pdf for correct rendering of all equations (an exact copy
of this file)

# Introduction

In this homework, you will be implementing ADMM algorithm to solve the
robust PCA problem. The idea behind robust PCA is as follows: given a
matrix *M* ∈ ℝ<sup>*n* × *p*</sup>, you want to decompose this matrix as
a sum of low-rank and sparse matrices. The standard PCA can be viewed as
a low-rank matrix decomposition of given *M* (no sparse matrix). In
contrast, the robust PCA uses an extra sparse matrix to allow for “few”
outliers - “few” elements of *M* that deviate from low rank structure
(hence extra sparse matrix, only “few” nonzero elements).

Given matrix *M*, Robust PCA is formulated as optimization problem
minimize<sub>*L*, *S*</sub>{∥*L*∥<sub>\*</sub>+*γ*∥*S*∥<sub>1</sub>}  s.t.  *L* + *S* = *M*.
Here ∥*L*∥<sub>\*</sub> is the nuclear norm of *L* (discussed in class),
∥*S*∥<sub>1</sub> = ∑<sub>*i*, *j*</sub>\|*s*<sub>*i*, *j*</sub>\| is
the sum of absolute values of the elements of matrix *S* (ℓ<sub>1</sub>
norm applied to vectorized *S*, i.e. lasso type penalty). The parameter
*γ* &gt; 0 controls the relative weights between the nuclear norm and
the sparsity penalties, and is chosen by the user.

The Robust PCA problem is already in the form required by ADMM, so the
implementation should consist of three updates: update of *L*, update of
*S* and update of dual variable *H* (matrix *η*, we consider scaled ADMM
for simplicity). You should discover that the first two updates
correspond to proximal operators of nuclear norm and ℓ<sub>1</sub> norm,
correspondingly, evaluated at specific points.

**Remark:** ADMM updates in vector form generalize easily to matrices by
substituting squared Frobenius norm ∥ ⋅ ∥<sub>*F*</sub><sup>2</sup>
instead of squared euclidean norm ∥ ⋅ ∥<sub>2</sub><sup>2</sup>, and the
matrix inner-product (trace operator tr(*A*<sup>⊤</sup>*B*)) instead of
the vector inner-product *a*<sup>⊤</sup>*b*.

# Derivation of ADMM updates

Your first task is to derive an explicit form of all three updates for
ADMM for Robust PCA problem above. To submit your derivations, please
**directly modify the README.Rmd** document in this folder (use R
markdown), and **knit** the document to update the corresponding pdf
accordingly.

**Augmented Lagrangian:** given *τ* &gt; 0 (*τ* = 1/*ρ*)

$$
L\_{\\tau}(L, S, \\nu) = \\{\|\|L\|\|\_{\*} + \\gamma \|\|S\|\|\_{1} + tr(V^{T}(L + S - M)) + \\frac{1}{2 \\tau} \|\|L + S - M\|\|^{2}\_{F} \\}
$$

Using the scaled version of ADMM updates, let *η* = *ν* ⋅ *τ*. Derive
explicit updates below.

**Update of L:**
$$
L^{t+1} = \\arg\\min\_{L}\\{\\mathcal{L}\_{\\tau}(L, S\_t, \\frac{\\eta\_t}{\\tau}) \\}, \\: (\\text{Where }\\nu\_t = \\frac{\\eta\_t}{\\tau})
$$
which is equivalent to (derive the explicit update and fill in below,
use more than one line if necessary)

**Update of S:**

$$
S^{t+1} = \\arg\\min\_{S} \\{\\mathcal{L}\_{\\tau}(L\_{t+1}, S, \\frac{\\eta\_t}{\\tau}) \\}, \\: (\\text{Where }\\nu\_t = \\frac{\\eta\_t}{\\tau})
$$
which is equivalent to (derive the explicit update and fill in below)

**Update of H:** (fill in, note that here *H* is a matrix version of
*η*)
$$
H^{t+1} = H^t + \\frac{1}{\\tau}(L\_{t+1} + S\_{t+1} - M)
$$

# Starter code

The starter code for four functions with detailed description is
provided in **ADMMfunctions.R**. The first three are helper functions,
whereas the last function is the full ADMM algorithm. Feel free to
create additional functions if needed, however please adhere to
specified format of the 4 provided functions.

-   There are two different soft-thresholding functions that you need to
    implement: 1) `soft` is a standard soft-thresholding but now
    generalized to work on matrix-valued input; 2) `soft_nuclear` is the
    soft-thresholding of singular values of a matrix that arises in the
    context of nuclear norm minimization.

The script **ApplyADMM.R** has a code to apply your function on a simple
synthetic example. The provided example is an illustration of how the
sparsity part allows to extract low-rank component even when the
original matrix is far from low-rank. You don’t need to add any
additional code to this script.

The script **TestADMM.R** is a place holder for your own tests. As
**ApplyADMM.R** is for exploratory purposes only, you still need to
extensively test your functions to make sure they are running correctly.

# Grading criteria

Your assignment will be graded based on

-   ADMM derivations *(10% of the grade)*

Filling out the ADMM updates in this .Rmd correctly.

-   correctness *(40% of the grade)*

Take advantage of objective function values over iterations as a way to
indirectly check the correctness of your function.

-   speed *(30% of the grade)*

There is not a lot of degrees of freedom in this assignment, but you
will get +5 bonus points if your **completely correct** code on 1st
submission is faster than mine (median your time/median mine time &lt;
0.9). My function in ApplyADMM.R script takes as median 37 sec to run on
my laptop for that specific example with that specific seed and
parameter selection.

-   code style/documentation *(10% of the grade)*

You need to comment different parts of the code so it’s clear what they
do, have good indentation, readable code with names that make sense. See
guidelines on R style, and posted grading rubric.

-   version control/commit practices *(10% of the grade)*

I expect you to start early on this assignment, and work gradually. You
want to commit often, have logically organized commits with short
description that makes sense. See guidelines on good commit practices,
and posted grading rubric.
