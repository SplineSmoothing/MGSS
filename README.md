# MGSS: A Multigrid Spline Smoothing Toolbox
Penalized splines is a popular method for function estimation under the assumption of smoothness and is well established for one- or two-dimensional covariates. The extension to multiple covariates is straightforward but suffers from exponentially increasing memory and computational complexity. This toolbox provides a matrix-free implementation of a conjugate gradient method (CG.R) for the regularized least squares problem resulting from tensor product B-spline smoothing with multivariate and scattered data as well as preconditioned version (MGCG.R) where a geometric multigrid preconditioner is applied. Further details can be found in [Siebenborn and Wagner, 2019](https://arxiv.org/abs/1901.00654).

## Manual
The manuals for the matrix-free CG-method (CG.R) and the matrix-free MGCG-method (MGCG.R) are provided.

### CG
After selecting the spline parameters, the transposed B-spline basis matrix and the curavture penalty were assembled for each spatial direction `p=1,...,P`:
```R
tPhi_list <- lapply(1:P, function(p) t( bspline_matrix(X[,p], m[p], q[p] ,Omega[[p]]) ) )     # spline matrices
Psi_list <- curvature_penalty(m, q, Omega)                                                    # curvature penalty
b <- MVP_khatrirao_rcpp(tPhi_list, y)                                                         # right-hand side vector
```
The coefficients of the spline basis functions are determined via the solution of a linear system which is achieved by the CG-method.
The key point is that the matrix-vector product with the coefficient matrix `A` are performed in a matrix-free manner, i.e. without explicitly assembling and storing the (too) large coefficient matrix but only the dimension specific matrices:
```R
Ad <- MVP_spline(tPhi_list, d) + lambda*MVP_penalty(Psi_list, d)
```

### MGCG
If the spatial dimension `P` is further increased, the CG-method will become computationally inefficient due to deteriorating condition of the system matrix. Therefore, an multigrid-like precondioner is implemented.
After selecting the spline parameters and the number of utilized grids, the transposed B-spline basis matrix and the curavture penalty were assembled for each spatial direction `p=1,...,P` and each grid level `g=1,...,G`:
```R
tPhi_list <- lapply(1:G, function(g) lapply(1:P, function(p) t( bspline_matrix(X[,p], m[[g]][p], q[p] ,Omega[[p]]) ) ) )    # spline matrices
Psi_list <- lapply(1:G, function(g)  curvature_penalty(m[[g]], q, Omega) )   # survature penalty
b <- MVP_khatrirao_rcpp(tPhi_list[[G]], y)      # right-hand side vector
```
The coefficients of the spline basis functions are determined via the solution of a linear system which is achieved by the preconditioned CG-method where the `v_cycle` is used as preconditioner:
```R
z <- v_cycle(tPhi_list, Psi_list, Rest, Prol, lambda, r, nu, w)     # apply MG as preconditioner
```
Note that all operations within the multigrid cycle, i.e. restriction, prolongation, and Jacobi-smoothing, are performed matrix-free such that none of the (too) large matrices will ever exist.

## Results
To evaluate the performance of the algorithms, they are tested on a rather simple data set in `P=3` dimensions
```R
P <- 3                                          # number of covariates
n <- 100000                                     # number of observations
X <- sapply(1:P, function(p) runif(n,0,1))      # covariates
t <- sapply( 1:n, function(i) -16*( (sum(X[i,]^2) / length(X[i,])) -0.5) )
fx <- 1 / ( 1 + exp(t) )
y <- fx + rnorm(n, 0, 0.1)                      # observations
```
Storing the full coefficient matrix of the underlying system requires approximately 30 GB of RAM, which is at the limit of the most computer systems.
The matrix-free approaches are, by construction, free of thoses memory limitations.
The results are as follows:

