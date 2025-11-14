//
// This Stan program defines an MCMC algorithm for estimating IID Fay Herriot model with variance smoothing
// 
//

functions {
  
   real pc_sigma_lpdf(real sigma, real U, real alpha){
    real lpdf;
    real theta;
    theta = -log(alpha)/U;
    lpdf = log(0.5*theta) + 3*log(sigma) - theta*sigma;
    return lpdf;
  }
}

data {
  int<lower=1> m;                  // number of areas
  vector[m] y;                     // direct estimates
  vector[m] v_hat;           // variance estimates
  
  //constants for Pearson approximation
  vector[m] Cons;
  vector[m] A;
  vector[m] B;
  vector[m] df;
  
  // Covariates for the mean model
  int<lower=1> p_mean;             // number of mean covariates
  matrix[m, p_mean] X;             // covariate matrix for all areas

  // Covariates for the variance model
  int<lower=1> p_var;              // number of variance covariates
  matrix[m, p_var] Z;              // covariate matrix for variance model
}

// The parameters accepted by the model. 
parameters {
  vector[p_mean] beta;             // coefficients for mean model
  vector[p_var] gamma;             // coefficients for variance model
  matrix[2, m] u_raw;   
  vector<lower=0>[2] sig_u;
  cholesky_factor_corr[2] L_u; // for correlation between random effects
}

transformed parameters {
  vector[m] theta;
  vector[m] log_sig2;
  vector<lower=0>[m] v_raw;
  matrix[2, m] u; 
  
  // Correlated random effects:
  u = (diag_pre_multiply(sig_u, L_u) * u_raw);


  // Mean and variance models
  theta = X * beta + u[1]';
  log_sig2 = Z * gamma + u[2]';
  
  v_raw = A.*v_hat./exp(log_sig2) + B;
}

// The model to be estimated.
model {
  target += pc_sigma_lpdf(sig_u[1] | 1, 0.01);
  target += pc_sigma_lpdf(sig_u[2] | 1, 0.01);
  L_u ~ lkj_corr_cholesky(2);        // weakly informative prior for correlation
  to_vector(u_raw) ~ normal(0, 1);   // standard normal base for REs
  
  beta ~ normal(0, 4);
  gamma[1] ~ normal(0, 4);   // intercept-like prior; adjust if Z includes intercept
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,4);
  }
 
  v_raw ~ chi_square(df);
  y ~ normal(theta,sqrt(exp(log_sig2).*Cons));
 
}

generated quantities {
  corr_matrix[2] Omega_u;
  real<lower=-1,upper=1> rho;
  
  Omega_u = multiply_lower_tri_self_transpose(L_u);  // implied correlation matrix
  rho = Omega_u[2,1];
}

