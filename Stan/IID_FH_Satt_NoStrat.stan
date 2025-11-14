//
// This Stan program 
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
  int<lower=1> m_data;                  // number of areas with data
  array[m_data] int data_areas;       // which area have data
  vector[m_data] y;                     // direct estimates
  vector[m_data] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  //constants for Satt approximation
  vector[m_data] Cons;
  vector[m_data] sum_q;
  vector[m_data] sum_q2;
  
  
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
//  cholesky_factor_corr[2] L_u; // for correlation between random effects
}

transformed parameters {
  vector[m] theta;
  vector[m_data] theta_data;
  vector[m] log_sig2;
  vector[m_data] log_sig2_data;
  vector<lower=0>[m_data] v_raw;
  vector<lower=0>[m_data] df;
  matrix[2, m] u; 
  
  // Correlated random effects:
 // u = (diag_pre_multiply(sig_u, L_u) * u_raw);
  u = diag_pre_multiply(sig_u, u_raw);

  // Mean and variance models
  theta = X * beta + u[1]';
  log_sig2 = Z * gamma + u[2]';
  
  for(i in 1:m_data){
    theta_data[i] = theta[data_areas[i]];
    log_sig2_data[i] = log_sig2[data_areas[i]];
  }
  
  v_raw = (sum_q./sum_q2).*v_hat_scaled./exp(log_sig2_data);
  df = (sum_q.*sum_q)./sum_q2;

}

// The model to be estimated.
model {
  target += pc_sigma_lpdf(sig_u[1] | 1, 0.01);
  target += pc_sigma_lpdf(sig_u[2] | 1, 0.01);
 // L_u ~ lkj_corr_cholesky(2);        // weakly informative prior for correlation
  to_vector(u_raw) ~ normal(0, 1);   // standard normal base for REs
  
  beta ~ normal(0, 4);
  gamma[1] ~ normal(0, 4);   // intercept-like prior; adjust if Z includes intercept
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,4);
  }
 
  v_raw ~ chi_square(df);
  y ~ normal(theta_data,sqrt(exp(log_sig2_data).*Cons));
 
}

//generated quantities {
//  corr_matrix[2] Omega_u;
//  real<lower=-1,upper=1> rho;
  
//  Omega_u = multiply_lower_tri_self_transpose(L_u);  // implied correlation matrix
//  rho = Omega_u[2,1];
//}

