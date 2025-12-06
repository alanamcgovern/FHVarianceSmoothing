
functions {
  
  real icar_lp(vector u, array[] int node1, array[] int node2) {
    // ICAR negative log density (intrinsic)
    int E = num_elements(node1);
    real out = 0;
    for (e in 1:E)
      out += (u[node1[e]] - u[node2[e]])^2;
    return -0.5 * out;
  }
  
   real pc_sigma_lp(real sigma, real U, real alpha){
    real lpdf;
    real theta;
    theta = -log(alpha)/U;
    lpdf = log(theta) - theta*sigma;
    return lpdf;
  }
}

data {
  int m;
  int m_data;
  array[m_data] int data_areas;
  vector[m_data] y;                     // direct estimates
  vector[m_data] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  int sig2_level;
  int len_sig2;
  array[m_data] int sig2_id;
  
  //constants for Satt approximation
  vector[m_data] df;
  
  // Covariates for the mean model
  int<lower=1> p_mean;             // number of mean covariates
  matrix[m, p_mean] X;             // covariate matrix for all areas

  // Covariates for the variance model
  int<lower=1> p_var;              // number of variance covariates
  matrix[len_sig2, p_var] Z;              // covariate matrix for variance model
  
  int<lower=0,upper=1> bym2_mean;   // 1 = BYM2 for mean model, 0 = IID
  int<lower=0,upper=1> bym2_var;    // 1 = BYM2 for variance model, 0 = IID
  
   // Adjacency (undirected) for ICAR admin1
  int<lower=0> N_edges1;            // number of edges
  array[N_edges1] int<lower=1, upper=m> node1_1;
  array[N_edges1] int<lower=1, upper=m> node1_2;
  real<lower=0> car_scale1;
  
   // Adjacency (undirected) for ICAR admin2
  int<lower=0> N_edges2;            // number of edges
  array[N_edges2] int<lower=1, upper=m> node2_1;
  array[N_edges2] int<lower=1, upper=m> node2_2;
  real<lower=0> car_scale2;
}

// The parameters accepted by the model. 
parameters {
  vector[p_mean] beta;             // coefficients for mean model
  vector[p_var] gamma;             // coefficients for variance model
  vector[m] u1;
  vector[len_sig2] u2;
  vector[m] s1_raw; // ICAR
  vector[len_sig2] s2_raw;
  vector<lower=0>[2] sig_u;
  real<lower=0,upper=1> phi1;
  real<lower=0,upper=1> phi2;
}

transformed parameters {
  vector[m] theta;
  vector[len_sig2] log_sig2;
  vector[m_data] theta_data;
  vector<lower=0>[m_data] sig2_data;
  vector<lower=0>[m_data] v_raw;
  vector[m] s1; // ICAR
  vector[len_sig2] s2;
  vector[m] b1;
  vector[len_sig2] b2;
  
  s1 = (s1_raw - mean(s1_raw))/sqrt(car_scale2);
  if(sig2_level ==2){
    s2 = (s2_raw - mean(s2_raw))/sqrt(car_scale2);
  }else{
    s2 = (s2_raw - mean(s2_raw))/sqrt(car_scale1);
  }
  
   // random effects:
  if(bym2_mean==1){
    b1 = sig_u[1]*(sqrt(phi1)*s1 + sqrt(1-phi1)*u1);
  }else{
    b1 = sig_u[1]*u1;
  }
  
  if(bym2_var==1){ // which scale to use for BYM2?
    b2 = sig_u[2]*(sqrt(phi2)*s2 + sqrt(1-phi2)*u2);
  }else{
    b2 = sig_u[2]*u2;
  }

  theta = X * beta + b1;
  log_sig2 = Z * gamma + b2;
  
  for(a in 1:m_data){
    theta_data[a] = theta[data_areas[a]];
    sig2_data[a] = exp(log_sig2[sig2_id[a]]);
  }
  
 v_raw = v_hat_scaled./sig2_data;
 
}

// The model to be estimated.
model {
  target += pc_sigma_lp(sig_u[1], 1, 0.01);
  target += pc_sigma_lp(sig_u[2], 1, 0.01);
  
  phi1 ~ beta(1,5);
  phi2 ~ beta(1,5);
  
  u1 ~ normal(0,1);
  u2 ~ normal(0,1);
  
  if(bym2_mean==1){
    target += icar_lp(s1_raw, node2_1,node2_2);
  }else{
    s1_raw ~ normal(0,1);
  }
  
  if(bym2_var==1 && sig2_level ==2){
    target += icar_lp(s2_raw, node2_1,node2_2);
  }else if(bym2_var==1 && sig2_level ==1){
    target += icar_lp(s2_raw, node1_1,node1_2);
  }else{
    s2_raw ~ normal(0,1);
  }
  
  beta[1] ~ normal(0, 2);
  if(p_mean>1){
     beta[2:p_mean] ~ normal(0,1);
  }
  gamma[1] ~ normal(-2, 2);
 // gamma[1] ~ normal(0, 2);
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,1);
  }
 
  v_raw ~ gamma(0.5*df,0.5);
  y ~ normal(theta_data,sqrt(sig2_data));

}
