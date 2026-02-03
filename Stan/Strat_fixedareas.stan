

data {
  int m;
  vector[m] y;                     // direct estimates
  vector[m] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
   vector[m] Cons;
  
  //constants for Satt approximation
  int lenq;
  array[lenq] int q_id;
  array[m] int q_start;
  array[m] int q_per_area;
  vector[lenq] q;
  vector[lenq] nu;
  
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
  vector[m] u1;
  vector[m] u2;
}

transformed parameters {
  vector[m] theta;
  vector[m] theta_drop;
  vector[m] log_sig2;
  vector<lower=0>[m] v;
  vector<lower=0>[m] v_raw;
  vector[lenq] delta;
  matrix[m,2] Q;
  vector[m] df;

  theta = X * beta + u1;
  log_sig2 = Z * gamma + u2;
  
  v = exp(log_sig2).*Cons;
  
  delta = square(nu*beta[2])./exp(log_sig2[q_id]);

  for(a in 1:m){
    Q[a,1] = sum(q[q_start[a]:(q_start[a]+q_per_area[a]-1)].*(1+delta[q_start[a]:(q_start[a]+q_per_area[a]-1)]));
    Q[a,2] = sum(square(q[q_start[a]:(q_start[a]+q_per_area[a]-1)]).*(1+2*delta[q_start[a]:(q_start[a]+q_per_area[a]-1)]));
  }
  
  df = square(Q[,1])./Q[,2];
  v_raw = (v_hat_scaled./v).*(Q[,1]./Q[,2]);

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
