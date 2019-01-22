

# matrix-free CG-method
my_CG = function(Phi_t_list, Psi_list, lambda, b, a=rep(0,length(b)), it_max=length(b), tol=10^(-6)){
  
  n_reg <- length(Psi_list)
  norm_b <- sqrt(sum(b^2)) 
  if(sum(a^2)!=0){
    base <- MVP_kr_Rcpp( Phi_t_list, MVP_krtrans_Rcpp(Phi_t_list,a) )
    pen <- Reduce( "+", lapply( 1:n_reg, function(p) MVP_kron_Rcpp(Psi_list[[p]],a) ) )
    Aa <- base + lambda*pen
    res <- b-Aa
  } else{
    res <- b
  }
  dir <- res
  norm_res <- sqrt(sum(res^2))
  for(k in 1:it_max){
    
    base <- MVP_kr_Rcpp( Phi_t_list, MVP_krtrans_Rcpp(Phi_t_list,dir) )
    pen <- Reduce("+" , lapply( 1:n_reg, function(p) MVP_kron_Rcpp(Psi_list[[p]],dir) ) )
    Adir <- base + lambda*pen
    alpha <- as.numeric(norm_res^2 / crossprod(dir,Adir))
    a <- a + alpha*dir
    res <- res - alpha*Adir
    norm_res_old <- norm_res
    norm_res <- sqrt(sum(res^2))
    if( (norm_res/norm_b) < tol){
      break
    }
    beta <- norm_res^2 / norm_res_old^2
    dir <- res + beta*dir
  }
  
  return(as.vector(a))
  
}


# matrix-free Jacobi
my_jacobi = function( Phi_t_list, Psi_d_list, diag_A, b, k_max, w=1, a=rep(0,length(b)), tol=10^(-8) ){
  
  # Parameter
  invD <- 1/diag_A
  W <- w*invD
  norm_b <- sqrt(sum(b^2))
  nreg <- length(Psi_d_list)
  
  # Jacobi Iteration
  if(sqrt(sum(a^2))!=0){
    base <- MVP_kr_Rcpp( Phi_t_list, MVP_krtrans_Rcpp(Phi_t_list,a) )
    pen <- Reduce( "+", lapply( 1:nreg, function(j) MVP_kron_Rcpp(Psi_d_list[[j]],a) ) )
    Aa <- base + pen
    res <- b-Aa
  } else{
    res <- b
  }
  for(k in 1:k_max){
    a <- a + W*res
    base <- MVP_kr_Rcpp( Phi_t_list, MVP_krtrans_Rcpp(Phi_t_list,a) )
    pen <- lapply( 1:nreg, function(j) MVP_kron_Rcpp(Psi_d_list[[j]],a) )
    Aa <- base + Reduce("+",pen)
    res <- b-Aa
    norm_res_rel <- sqrt(sum(res^2))/norm_b
    if(norm_res_rel < tol){
      break
    }
  }
  return(a)

}


# prolongation matrix I_{g-1}^{g} (coarse to fine)
my_prolongation = function(K_coarse, K_fine, q){
  
  g <- q+1
  values <- (1/(2^q))*choose(g, 0:g)   # Pascals Triangle
  k_zeros <- K_fine+2*q-length(values)
  S <- sapply(1:K_coarse, function(j)  c( rep(0,2*(j-1)) , values , rep(0,k_zeros-(2*(j-1))) ) )
  P <- S[g:(K_fine+q),]
  return(P)
  
}


# restriction matrix I_{g}^{g-1} (fine to coarse)
my_restriction = function(K_fine, K_coarse, q){
  
  R <- t( my_prolongation(K_coarse, K_fine, q) )
  return(R)
  
}


# matrix-free v-cycle (with Jacobi smoother and CG coarse grid solver)
my_vcycle = function( Phi_t_list, Psi_list, diags, Rest, Prol, b, nu, w=1, a=rep(0,length(b)) ){
  
  g <- length(Phi_t_list)
  nreg <- length(Psi_list[[g]])
  
  if(g==1){
    
    z <- my_CG(Phi_t_list[[g]], Psi_list[[g]], lambda=1, b)
    #z <- forwardsolve(t(U_chol), as.vector(b) )
    #return( backsolve(U_chol, z) )
    
  } else{
    
    a <- my_jacobi(Phi_t_list[[g]], Psi_list[[g]], diags[[g]], b, nu[1], w, a )
    base <- MVP_kr_Rcpp( Phi_t_list[[g]], MVP_krtrans_Rcpp(Phi_t_list[[g]],a) )
    pen <- lapply( 1:nreg, function(j) MVP_kron_Rcpp(Psi_list[[g]][[j]],a) )
    Aa <- base +  Reduce( "+",pen)
    e <- b - Aa
    r <- MVP_kron_Rcpp(Rest[[g-1]], e)
    e <- my_vcycle( Phi_t_list[1:(g-1)], Psi_list[1:(g-1)], diags[1:(g-1)], Rest[1:(g-1)], Prol[1:(g-1)], r, nu, w )
    a <- a + MVP_kron_Rcpp( Prol[[g-1]], e )
    a <- my_jacobi(Phi_t_list[[g]], Psi_list[[g]], diags[[g]], b, nu[2], w, a )
    
  }
  
  return(as.vector(a))
  
}

