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
  
  //constants for Satt approximation
  int len_obs;
  vector[len_obs] urban_id;
  vector[len_obs] n_vec;
  vector[len_obs] wt_vec;
  matrix[m_data,2] Cons;
  array[m_data] int start_id;
  array[m_data] int obs_per_area;
  array[m_data] int ranks;
  int ncol_C;
  matrix[len_obs,ncol_C] C_long;
  
  //real kappa;
 
  // Covariates for the mean model
  int<lower=1> p_mean;             // number of mean covariates
  matrix[m, p_mean] X;             // covariate matrix for all areas

  // Covariates for the variance model
  int<lower=1> p_var;              // number of variance covariates
  matrix[m, p_var] Z;              // covariate matrix for variance model
  
  int<lower=0,upper=1> bym2_mean;   // 1 = BYM2 for mean model, 0 = IID
  int<lower=0,upper=1> bym2_var;    // 1 = BYM2 for variance model, 0 = IID
  
   // Adjacency (undirected) for ICAR admin2
  int<lower=0> N_edges;            // number of edges
  array[N_edges] int<lower=1, upper=m> node_1;
  array[N_edges] int<lower=1, upper=m> node_2;
  real<lower=0> car_scale;
}

// The parameters accepted by the model. 
parameters {
  vector[p_mean] beta;             // coefficients for mean model
  vector[p_var] gamma;             // coefficients for variance model
  vector[m] u1;
  vector[m] u2;
  vector[m] s1_raw; // ICAR
  vector[m] s2_raw;
  vector<lower=0>[2] sig_u;
  real<lower=0,upper=1> phi1;
  real<lower=0,upper=1> phi2;
  real<lower=-1> kappa;
}

transformed parameters {
  
  real kappa_trans = log(kappa+1);
  
  vector[m] s1 = (s1_raw - mean(s1_raw))/sqrt(car_scale);
  vector[m] s2 = (s2_raw - mean(s2_raw))/sqrt(car_scale);
  
   // random effects:
  vector[m] b1;
  vector[m] b2;
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

  vector[m] theta = X * beta + b1;
  vector[m] log_v = Z * gamma + b2;
  
  vector[m_data] theta_data;
  vector<lower=0>[m_data] v_data;
  matrix[m_data,2] sums;
  
  for(a in 1:m_data){
    vector[obs_per_area[a]] wt_tmp= wt_vec[start_id[a]:(start_id[a] + obs_per_area[a]-1)];
    vector[obs_per_area[a]] n_tmp= n_vec[start_id[a]:(start_id[a] + obs_per_area[a]-1)];
    vector[obs_per_area[a]] urban_tmp= urban_id[start_id[a]:(start_id[a] + obs_per_area[a]-1)];
    vector[obs_per_area[a]] S_sqrt_vec= sqrt((1+kappa*urban_tmp)./n_tmp);
    matrix[obs_per_area[a],obs_per_area[a]] C= C_long[start_id[a]:(start_id[a] + obs_per_area[a]-1),1:obs_per_area[a]];
    
    tuple(matrix[obs_per_area[a],obs_per_area[a]], vector[obs_per_area[a]]) ev = eigendecompose_sym(diag_matrix(S_sqrt_vec)*C*diag_matrix(S_sqrt_vec));
    matrix[obs_per_area[a],ranks[a]] U=ev.1[,(obs_per_area[a]-ranks[a]+1):obs_per_area[a]];
    vector[ranks[a]] q=ev.2[(obs_per_area[a]-ranks[a]+1):obs_per_area[a]];
    
    vector[ranks[a]] delta = sum(square(wt_tmp).*n_tmp)/square(sum(wt_tmp.*n_tmp))*square(U'*diag_matrix(1/S_sqrt_vec)*urban_tmp*beta[2])./exp(log_v[data_areas[a]]);

    sums[a,1] = sum(q.*(1+delta));
    sums[a,2] = sum(square(q).*(1+2*delta));
    
    theta_data[a] = theta[data_areas[a]];
    v_data[a] = exp(log_v[data_areas[a]]);
  }
  
  vector[m_data] df = square(sums[,1])./sums[,2];
  vector<lower=0>[m_data] v_raw = v_hat_scaled./v_data.*sums[,1]./sums[,2].*(Cons[,1]+kappa*Cons[,2]);

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
    target += icar_lp(s1_raw, node_1,node_2);
  }else{
    s1_raw ~ normal(0,1);
  }
  
  if(bym2_var==1){
    target += icar_lp(s2_raw, node_1,node_2);
  }else{
    s2_raw ~ normal(0,1);
  }
  
  kappa_trans ~ normal(0,0.5);
  
  beta ~ normal(0, 2);
  
  gamma[1] ~ normal(-3, 2);
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,1);
  }
 
  v_raw ~ gamma(0.5*df,0.5);
  y ~ normal(theta_data,sqrt(v_data));

}
