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
  
   real icar_normal_lpdf(vector u,int N,array[] int node1,array[] int node2){ 
    return -0.5 * dot_self(u[node1]-u[node2]) +normal_lpdf(sum(u)|0,0.001 * N); 
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
  vector[m_data] df;
  
  // Covariates for the mean model
  int<lower=1> p_mean;             // number of mean covariates
  matrix[m, p_mean] X;             // covariate matrix for all areas

  // Covariates for the variance model
  int<lower=1> p_var;              // number of variance covariates
  matrix[m, p_var] Z;              // covariate matrix for variance model
  
   // Adjacency (undirected) for ICAR
  int<lower=0> N_edges;            // number of edges
  array[N_edges] int<lower=1, upper=m> node1;
  array[N_edges] int<lower=1, upper=m> node2;

  // Scaling factor to put ICAR on the BYM2 scale
  // (typical sd of raw ICAR field; see Riebler et al. 2016)
  real<lower=0> car_scale;
}

// The parameters accepted by the model. 
parameters {
  vector[p_mean] beta;             // coefficients for mean model
  vector[p_var] gamma;             // coefficients for variance model
   vector[m] u1; // IID
  vector[m] u2;
  vector[m] s1; // ICAR
 // vector[m] s2;
  vector<lower=0>[2] sig_u;
  //vector<lower=0>[2] phi;
  real<lower=0,upper=1> phi;
}

transformed parameters {
  vector[m] theta;
  vector[m_data] theta_data;
  vector[m] log_sig2;
  vector[m_data] log_sig2_data;
   vector[m] b1;
  vector[m] b2;
  
  //  random effects:
   b1 = sig_u[1]*(sqrt(phi/car_scale)*s1 + sqrt(1-phi)*u1);
  b2 = sig_u[2]*u2;
  
  //b1 = sig_u[1]*(sqrt(phi[1]/car_scale)*s1 + sqrt(1-phi[1])*u1);
  //b2 = sig_u[2]*(sqrt(phi[2]/car_scale)*s2 + sqrt(1-phi[2])*u2);

 
  // Mean and variance models
  theta = X * beta + b1;
  log_sig2 = Z * gamma + b2;
  
  for(i in 1:m_data){
    theta_data[i] = theta[data_areas[i]];
    log_sig2_data[i] = log_sig2[data_areas[i]];
  }

}

// The model to be estimated.
model {
  target += pc_sigma_lpdf(sig_u[1] | 1, 0.01);
  target += pc_sigma_lpdf(sig_u[2] | 1, 0.01);
 phi ~ beta(0.5,0.5);
  
 // phi[1] ~ beta(0.5,0.5);
  //phi[2] ~ beta(0.5,0.5);
  
  u1 ~ normal(0, 1);
  u2 ~ normal(0, 1);
  
  target += icar_normal_lpdf(s1 | m,node1,node2);
  //target += icar_normal_lpdf(s2 | m,node1,node2);
  
  beta ~ normal(0, 2);
  gamma[1] ~ normal(0, 2);   // intercept-like prior; adjust if Z includes intercept
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,2);
  }
 
  v_hat_scaled ~ gamma(0.5*df,0.5);
  y ~ normal(theta_data,sqrt(exp(log_sig2_data).*Cons));
 
}
