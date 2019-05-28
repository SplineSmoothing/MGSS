

#####----------------------------------------------------------------------------------------------------------------
#####  prolongation matrix I_{g-1}^{g} (coarse to fine)
prolongation_matrix = function(K_coarse, K_fine, q){
  
  g <- q+1
  values <- (1/(2^q))*choose(g, 0:g)   # Pascals Triangle
  k_zeros <- K_fine+2*q-length(values)
  S <- sapply(1:K_coarse, function(j)  c( rep(0,2*(j-1)) , values , rep(0,k_zeros-(2*(j-1))) ) )
  P <- S[g:(K_fine+q),]
  return(P)
  
}


#####----------------------------------------------------------------------------------------------------------------
#####  restriction matrix I_{g}^{g-1} (fine to coarse)
restriction_matrix = function(K_fine, K_coarse, q){
  
  R <- t( prolongation_matrix(K_coarse, K_fine, q) )
  return(R)
  
}


#####----------------------------------------------------------------------------------------------------------------
#####  matrix-free CG method
solve_CG = function(tPhi_list, Psi_list, lambda, b, alpha=rep(0,length(b)), max_iter=length(b), tol=10^(-6)){
  
  norm_b <- sqrt(sum(b^2))
  if( max(alpha^2)!=0 ){
    Aalpha <- MVP_spline(tPhi_list, alpha) + lambda*MVP_penalty(Psi_list, alpha)
    r <- b-Aalpha
  } else{
    r <- b
  }
  d <- r
  r_square <- crossprod(r)
  for(i in 1:max_iter){
    
    Ad <- MVP_spline(tPhi_list, d) + lambda*MVP_penalty(Psi_list, d)
    t <- as.numeric( r_square / crossprod(d, Ad) )
    alpha <- alpha+t*d
    r_new <- r-t*Ad
    r_new_square <- crossprod(r_new)
    if(sqrt(r_new_square) < tol){   
      break
    }
    beta <- as.numeric( r_new_square / r_square )   
    r <- r_new                                      
    r_square <- r_new_square                        
    d <- r+beta*d                 
    
  }
  
  return(as.vector(alpha))
  
}


#####----------------------------------------------------------------------------------------------------------------
#####  matrix-free Jacobi as smoothing iteration
smooth_jacobi = function(tPhi_list, Psi_list, lambda, b, k_max, w=1, a=rep(0,length(b)), tol=10^(-8) ){
  
  # Parameter
  diag_A <- diag_khatrirao_rcpp(tPhi_list) + lambda*rowSums(sapply( 1:length(Psi_list), function(j) diag_kronecker_rcpp(Psi_list[[j]]) ) )
  invD <- 1/diag_A
  W <- w*invD
  norm_b <- sqrt(sum(b^2))
  nreg <- length(Psi_list)
  
  # Jacobi Iteration
  if(sqrt(sum(a^2))!=0){
    Aa <- MVP_spline(tPhi_list, a) + lambda*MVP_penalty(Psi_list, a)
    res <- b-Aa
  } else{
    res <- b
  }
  for(k in 1:k_max){
    a <- a + W*res
    Aa <- MVP_spline(tPhi_list, a) + lambda*MVP_penalty(Psi_list, a)
    res <- b-Aa
    if(sqrt(sum(res^2))/norm_b < tol){
      break
    }
  }
  return(a)

}


#####----------------------------------------------------------------------------------------------------------------
#####  matrix-free v-cycle (with Jacobi smoother and CG coarse grid solver)
v_cycle = function(tPhi_list, Psi_list, Rest, Prol, lambda, b, nu, w=1, a=rep(0,length(b)) ){
  
  g <- length(tPhi_list)
    if(g==1){
    
    z <- solve_CG(tPhi_list[[g]], Psi_list[[g]], lambda, b)

  } else{
    
    a <- smooth_jacobi(tPhi_list[[g]], Psi_list[[g]], lambda, b, nu[1], w, a )
    Aa <- MVP_spline(tPhi_list[[g]], a) + lambda*MVP_penalty(Psi_list[[g]], a)
    e <- b - Aa
    r <- MVP_kronecker_rcpp(Rest[[g-1]], e)
    e <- v_cycle(tPhi_list[1:(g-1)], Psi_list[1:(g-1)], Rest[1:(g-1)], Prol[1:(g-1)], lambda, r, nu, w )
    a <- a + MVP_kronecker_rcpp( Prol[[g-1]], e )
    a <- smooth_jacobi(tPhi_list[[g]], Psi_list[[g]], lambda, b, nu[2], w, a )
    
  }
  
  return(as.vector(a))
  
}

