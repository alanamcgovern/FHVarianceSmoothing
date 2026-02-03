

data {
  int m;
  vector[m] y;                     // direct estimates
  vector[m] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  //constants 
  vector[m] df;
  vector[m] Cons;
  
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
  vector[m] u1; // fixed area effect for mean
  vector[m] u2; // fixed area effect for variance
}

transformed parameters {
  vector[m] theta;
  vector[m] log_sig2;
  vector<lower=0>[m] v;
  vector<lower=0>[m] v_raw;
  
  theta = X * beta + u1;
  log_sig2 = Z * gamma + u2;
  
  v = exp(log_sig2).*Cons;
  v_raw = v_hat_scaled./v;
 
}

// The model to be estimated.
model {
  
  u1 ~ normal(0,1);
  u2 ~ normal(0,1);
  
  beta ~ normal(0, 2);
  
  gamma[1] ~ normal(0.5, 0.5);
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,1);
  }
 
  v_raw ~ gamma(0.5*df,0.5);
  y ~ normal(theta,sqrt(v));

}
