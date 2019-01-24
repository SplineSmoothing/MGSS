# MGSS - A Multigrid Spline Smoothing Toolbox
This package provides two matrix-free algorithms to determine smoothing splines with multiple covariates.
Theses algorithms are:  
* Matrix-free CG-method (CG.R)
* Matrix-free MGCG-method (MGCG.R)

To test both of the algorithms in multiple dimensions, we provide test data in P dimensions (P=2,3,4) that is generated from a disturbed sigmoid function.
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

### For MGCG.R
G <- 5    # number of grids
m <- lapply( 1:G, function(g) rep(2^g-1,P) )
```
and assemble the required matrices and vectors:
```R
### For CG.R:
Phi_t_list <- lapply(1:P, function(p) t(my_bs_matrix(X[,p],m[p],q[p],Omega[[p]])) )
Psi_list <- my_TP_regularization(m,q,Omega)
b <- MVP_kr_Rcpp(Phi_t_list,y)

### For MGCG.R
Phi_t_list <- lapply( 1:G, function(g) lapply( 1:P, function(p) t(my_bs_matrix(X[,p],m[[g]][p],q[p],Omega[[p]])) ) )
Psi_list <- lapply( 1:G, function(g) my_TP_regularization(m[[g]],q,Omega) )
b <- MVP_kr_Rcpp(Phi_t_list[[G]],y)
```
Note that for the CG.R method, other basis and penalty matrices can be used as well.
We manually set the smoothing parameter:
```R
lambda <- 0.2
```
and finally solve the considered large-scale linear system with the matrix-free CG algorithm:
```R
alpha <- my_CG(Phi_t_list, Psi_list, lambda, b)
```
or the the repetetive aplication of the matrix-free v-cycle:
```R
alpha <- my_vcycle( Phi_t_list, Psi_list, diags, Rest, Prol, b, nu, w, alpha)
```

