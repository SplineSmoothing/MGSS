
###
library(Rcpp)
library(gaussquad)
library(combinat)
source("../functions/RLS_functions.R")
source("../functions/MG_functions.R")
sourceCpp("../functions/Rcpp_functions.cpp")

### data
P <- 2
data <- read.csv( paste("../data/data_",P,"D.csv",sep=""), dec="," )
X <- data[,1:P]
y <- data$y
n <- length(y)

### parameter
m <- rep(36,P)
q <- rep(3,P)
#J <- m+q+1
#K <- prod(J)
Omega <- lapply( 1:P, function(p) c(min(X[,p]),max(X[,p])) )

### basis matrices Phi_p^T
Phi_t_list <- lapply(1:P, function(p) t(my_bs_matrix(X[,p],m[p],q[p],Omega[[p]])) )

### thin plate penalty matrices Psi_rp^p
Psi_list <- my_TP_regularization(m,q,Omega)

### right-hand side vector
b <- MVP_kr_Rcpp(Phi_t_list,y)

### smoothing parameter (manually)
lambda <- 0.2

### CG method
alpha <- my_CG(Phi_t_list, Psi_list, lambda, b)


