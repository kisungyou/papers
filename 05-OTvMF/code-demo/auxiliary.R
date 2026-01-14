# auxiliary functions to be called

# aux_density2     : generate a 2-dimensional density line graph for vMF on S^1
# aux_density2_mix   - similar job but a bit more
# aux_quickbary    : given two vMF, interpolate quickly
# aux_pdist_L2     : Monte Carlo distance calculation
# aux_movMF_exact  : fit the movMF data exactly and return information
# aux_movMF_reduce : fit a big model, reduce it, and return information
# aux_vMF_estimate : estimate two parameters

library(movMF)
library(maotai)

# estimate vMF parameters -------------------------------------------------
aux_vMF_estimate <- function(X){
  fitted = movMF::movMF(X, k=1)
  vec_theta = as.vector(fitted$theta)
  kappa = sqrt(sum(vec_theta^2))
  mu = vec_theta/kappa
  
  output = list(mu=mu, kappa=kappa)
  return(output)
}

# aux_movMF_reduce --------------------------------------------------------
aux_movMF_reduce <- function(X, big.K=10, small.K=2, use.greedy = TRUE, use.hclust=TRUE){
  # fit the model
  fit_big_obj <- movMF::movMF(X, k=big.K)
  
  # convert into maotai-compatible format and extract
  fit_big_maotai <- maotai::movMF_convert(fit_big_obj)
  big_mu <- fit_big_maotai$means
  big_kappa <- fit_big_maotai$concentrations
  big_weight <- fit_big_maotai$weights
  
  # reduce the model
  if (use.greedy){
    fit_small_maotai <- maotai::movMF_reduce_greedy(big_mu, 
                                                    big_kappa,
                                                    big_weight, 
                                                    small.K) 
  } else {
    if (use.hclust){
      fit_small_maotai <- maotai::movMF_reduce_partitional(big_mu, 
                                                         big_kappa,
                                                         big_weight, 
                                                         small.K, method="hclust")
    } else {
      fit_small_maotai <- maotai::movMF_reduce_partitional(big_mu, 
                                                         big_kappa,
                                                         big_weight, 
                                                         small.K, method="kmedoids")
    }
  }
  
  # extract small model pars
  small_mu = fit_small_maotai$means
  small_kappa = fit_small_maotai$concentrations
  small_weight = fit_small_maotai$weights
  
  # evaluate the information
  fit_info <- maotai::movMF_info(X, small_mu, small_kappa, small_weight)
  return(fit_info)
}

# aux_movMF_exact ---------------------------------------------------------
aux_movMF_exact <- function(X, K=1){
  # fit the model
  fit_obj <- movMF::movMF(X, k=K)
  
  # convert into maotai-compatible format and extract
  fit_maotai <- maotai::movMF_convert(fit_obj)
  fit_mu = fit_maotai$means
  fit_kappa = as.vector(fit_maotai$concentrations)
  fit_weight = fit_maotai$weights
  
  # evaluate the information
  fit_info <- maotai::movMF_info(X, fit_mu, fit_kappa, fit_weight)
  
  return(fit_info)
}


# aux_crossdist_WL --------------------------------------------------------
aux_crossdist_WL <- function(mu1, kappa1, mu2, kappa2, scaler=1){
  M = length(kappa1)
  N = length(kappa2)
  
  output = array(0,c(M,N))
  for (m in 1:M){
    for (n in 1:N){
      dist1 <- base::sum((as.vector(mu1[m,]) - as.vector(mu2[n,]))^2)
      dist2 <- (base::ncol(mu1)-1)*scaler*((1/sqrt(kappa1[m]) - 1/sqrt(kappa2[n])))^2
      output[m,n] = sqrt(dist1+dist2)
    }
  }
  return(output)
}

# aux_pdist_WL ------------------------------------------------------------
aux_pdist_WL <- function(means, concentrations, scaler=1){
  n = length(concentrations)
  d = ncol(means)
  
  output = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      dist1 <- base::sum((as.vector(means[i,])-as.vector(means[j,]))^2)
      #dist1 <- acos(sum(as.vector(means[i,])*as.vector(means[j,])))^2
      dist2 <- (base::ncol(means)-1)*scaler*((1/sqrt(concentrations[i]) - 1/sqrt(concentrations[j])))^2
      
      output[i,j] <- output[j,i] <- sqrt(dist1+dist2)
    }
  }
  return(output)
}


# aux_crossdist_L2 --------------------------------------------------------
aux_crossdist_L2 <- function(mu1, kappa1, mu2, kappa2){
  d = base::ncol(mu1)
  
  # Uniform sampling on the unit sphere
  sample_sphere <- function(n, d) {
    mat <- matrix(rnorm(n * d), nrow = n, ncol = d)
    mat <- mat / sqrt(rowSums(mat^2))  # Normalize each row to have unit norm
    return(mat)
  }
  
  # Compute the vMF density
  vmf_density <- function(x, mu, kappa) {
    C_d <- (kappa^(d/2 - 1)) / ((2 * pi)^(d/2) * besselI(kappa, d/2 - 1, expon.scaled = TRUE))
    exp(kappa * (x %*% mu)) * C_d
  }
  
  vmf_l2_distance <- function(mu1, kappa1, mu2, kappa2, n_samples = 100) {
    d <- length(mu1)  # Dimensionality
    
    # Uniformly sample from the unit sphere
    X <- sample_sphere(n_samples, d)
    
    # Compute densities
    p1 <- apply(X, 1, function(x) vmf_density(x, mu1, kappa1))
    p2 <- apply(X, 1, function(x) vmf_density(x, mu2, kappa2))
    
    # Compute volume of the sphere
    sphere_volume <- 2 * pi^(d/2) / gamma(d/2)
    
    # Monte Carlo estimate with volume correction
    integral_estimate <- mean((p1 - p2)^2) * sphere_volume
    
    return(sqrt(integral_estimate))  # L2 distance
  }
  
  M = length(kappa1)
  N = length(kappa2)
  output = array(0,c(M,N))
  for (m in seq_len(M)){
    for (n in seq_len(N)){
      output[m,n] = vmf_l2_distance(mu1[m,], kappa1[m], mu2[n,], kappa2[n])
    }
  }
  return(output)
}

# aux_pdist_L2 ------------------------------------------------------------
aux_pdist_L2 <- function(means, concentrations){
  d <- base::ncol(means)
  
  # Uniform sampling on the unit sphere
  sample_sphere <- function(n, d) {
    mat <- matrix(rnorm(n * d), nrow = n, ncol = d)
    mat <- mat / sqrt(rowSums(mat^2))  # Normalize each row to have unit norm
    return(mat)
  }
  
  # Compute the vMF density
  vmf_density <- function(x, mu, kappa) {
    C_d <- (kappa^(d/2 - 1)) / ((2 * pi)^(d/2) * besselI(kappa, d/2 - 1, expon.scaled = TRUE))
    exp(kappa * (x %*% mu)) * C_d
  }
  
  vmf_l2_distance <- function(mu1, kappa1, mu2, kappa2, n_samples = 100) {
    d <- length(mu1)  # Dimensionality
    
    # Uniformly sample from the unit sphere
    X <- sample_sphere(n_samples, d)
    
    # Compute densities
    p1 <- apply(X, 1, function(x) vmf_density(x, mu1, kappa1))
    p2 <- apply(X, 1, function(x) vmf_density(x, mu2, kappa2))
    
    # Compute volume of the sphere
    sphere_volume <- 2 * pi^(d/2) / gamma(d/2)
    
    # Monte Carlo estimate with volume correction
    integral_estimate <- mean((p1 - p2)^2) * sphere_volume
    
    return(sqrt(integral_estimate))  # L2 distance
  }
  
  N = length(concentrations)
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      output[i,j] = vmf_l2_distance(means[i,], concentrations[i], means[j,], concentrations[j])
      output[j,i] <- output[i,j]
    }
  }
  
  return(output)
}


# aux_quickbary -----------------------------------------------------------
aux_quickbary <- function(mu1, kappa1, mu2, kappa2, alpha){
  out_mu = (1-alpha)*mu1 + alpha*mu2
  out_mu = out_mu/sqrt(sum(out_mu^2))
  
  out_kappa = ((1-alpha)/sqrt(kappa1) + (alpha/sqrt(kappa2)))^(-2)
  
  return(list(mu=out_mu,kappa=out_kappa))
}

# aux_density2 ------------------------------------------------------------
#   theta: angle values in a vector format
#   mu: mean direction
#   kappa: concentration
aux_density <- function(theta, mu, kappa){
  # parameter check
  par_mu = mu/sqrt(sum(mu^2))
  par_kappa = as.double(kappa)
  
  # convert angles to unit vectors on S^1
  X <- cbind(cos(theta), sin(theta))
  
  # compute vMF density
  c_1 <- par_kappa / (2 * pi * besselI(par_kappa, 0, expon.scaled = TRUE)) * exp(par_kappa)
  density <- rep(0, length(theta))
  for (i in 1:length(theta)){
    density[i] = c_1*exp(par_kappa*sum(par_mu*as.vector(X[i,])))
  }
  return(density)
}

# assume equal weights
aux_density2_mix <- function(theta, locations, concentrations){
  K = length(concentrations)
  mat_densities = array(0,c(length(theta),K))
  
  # convert angles to unit vectors on S^1
  X <- cbind(cos(theta), sin(theta))
  
  # for each vMF, compute the densities
  for (k in 1:K){
    # parameter check
    par_mu = as.vector(locations[k,]); par_mu = par_mu/sqrt(sum(par_mu^2))
    par_kappa = as.double(concentrations[k])
    
    # compute
    c_1 <- par_kappa / (2 * pi * besselI(par_kappa, 0, expon.scaled = TRUE)) * exp(par_kappa)
    for (n in 1:length(theta)){
      mat_densities[n,k] = c_1*exp(par_kappa*sum(par_mu*as.vector(X[n,])))
    }
  }
  
  # weighted sum and normalization
  density_vec  <- rowSums(mat_densities)/4
  density_norm <- density_vec/max(density_vec)
  
  # convert angles to cartesian coordinates
  x_coords <- cos(theta)*(1 + 0.3*density_norm)
  y_coords <- sin(theta)*(1 + 0.3*density_norm)
  
  # return the output
  output = list()
  output$coords = cbind(x_coords, y_coords)
  #output$whichK = apply(mat_densities, 1, which.max)
  return(output)
}

aux_density2 <- function(theta, mu, kappa, scaler=1){
  # parameter check
  par_mu = mu/sqrt(sum(mu^2))
  par_kappa = as.double(kappa)
  
  # convert angles to unit vectors on S^1
  X <- cbind(cos(theta), sin(theta))
  
  # compute vMF density
  c_1 <- par_kappa / (2 * pi * besselI(par_kappa, 0, expon.scaled = TRUE)) * exp(par_kappa)
  density <- rep(0, length(theta))
  for (i in 1:length(theta)){
    density[i] = c_1*exp(par_kappa*sum(par_mu*as.vector(X[i,])))
  }
  
  # normalize density
  density_normalized <- density/base::max(density)

  # convert angles to cartesian coordinates
  x_coords <- cos(theta)*(1 + 0.3*density_normalized*scaler)
  y_coords <- sin(theta)*(1 + 0.3*density_normalized*scaler)

  # return the output
  return(cbind(x_coords, y_coords))
}


