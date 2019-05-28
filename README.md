# MGSS: A Multigrid Spline Smoothing Toolbox
Penalized splines is a popular method for function estimation under the assumption of smoothness and is well established for one- or two-dimensional covariates. The extension to multiple covariates is straightforward but suffers from exponentially increasing memory and computational complexity. This toolbox provides a matrix-free implementation of a geometric multigrid preconditioned conjugate gradient method for the regularized least squares problem resulting from tensor product B-spline smoothing with multivariate and scattered data. Further details can be found in 




This package provides matrix-free algorithms to determine smoothing splines with multiple covariates.
The algorithms are:  
* Matrix-free CG-method (CG.R)
* Matrix-free MG-method (MG.R)
* Matrix-free MGCG-method (MGCG.R)

To test the algorithms in multiple dimensions, we provide test data in P dimensions (P=2,3,4) that is generated from a disturbed sigmoid function.
Before running one of the algorithms install the required packages:
```R
install.packages(c("Rcpp","gaussquad","combinat","polynom","orthopolynom"))
```
and set the dimesion:
```R
P <- 3 #or 2 or 4
```
As a next step, fix the spline parameters:
```R
q <- rep(3,P)   # spline degree in spatial direction

### For CG.R
m <- rep(36,P)  # number of knots in spatial direction

### For MG.R and MGCG.R
G <- 5    # number of grids
m <- lapply( 1:G, function(g) rep(2^g-1,P) )
```
and assemble the required matrices and vectors:
```R
### For CG.R:
Phi_t_list <- lapply(1:P, function(p) t(my_bs_matrix(X[,p],m[p],q[p],Omega[[p]])) )
Psi_list <- my_TP_regularization(m,q,Omega)
b <- MVP_kr_Rcpp(Phi_t_list,y)

### For MG.R and MGCG.R
Phi_t_list <- lapply( 1:G, function(g) lapply( 1:P, function(p) t(my_bs_matrix(X[,p],m[[g]][p],q[p],Omega[[p]])) ) )
Psi_list <- lapply( 1:G, function(g) my_TP_regularization(m[[g]],q,Omega) )
b <- MVP_kr_Rcpp(Phi_t_list[[G]],y)
```
Note that for the CG.R method, other basis and penalty matrices can be used as well.
We manually set the smoothing parameter:
```R
lambda <- 0.2
```
and finally solve the considered large-scale linear system with one of the preoposed algorithms.
