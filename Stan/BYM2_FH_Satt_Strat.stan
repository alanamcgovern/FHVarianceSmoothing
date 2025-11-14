//
// This Stan program 
// 
//

functions {
  
   real icar_normal_lpdf(vector u,int N,array[] int node1,array[] int node2){ 
    return -0.5 * dot_self(u[node1]-u[node2]) +normal_lpdf(sum(u)|0,0.001 * N); 
  }
  
   real pc_sigma_lpdf(real sigma, real U, real alpha){
    real lpdf;
    real theta;
    theta = -log(alpha)/U;
    lpdf = log(0.5*theta) + 3*log(sigma) - theta*sigma;
    return lpdf;
  }
  
  vector interpolate_vec(real kappa, vector kappa_seq, matrix vals) {
    int K = num_elements(kappa_seq);
    vector[K-1] is_lower;
    vector[K-1] w;

    for (j in 1:(K-1)) {
        if (kappa >= kappa_seq[j] && kappa <= kappa_seq[j+1]) {
            real weight = (kappa - kappa_seq[j]) /
                          (kappa_seq[j+1] - kappa_seq[j]);
            return weight * vals[,j+1] + (1 - weight) * vals[,j];
        }
    }
    // if outside range, clamp to endpoints
    if (kappa < kappa_seq[1]) return vals[,1];
    if (kappa > kappa_seq[K]) return vals[,K];
    else return vals[,1];

}

}

data {
  int<lower=1> m;                  // number of areas
  int<lower=1> m_data;                  // number of areas with data
  array[m_data] int data_areas;       // which area have data
  
  int<lower=m_data> len_q;
  int len_kappa_seq;
  vector[len_kappa_seq] kappa_seq;
  array[m_data] int q_per_area; 
  array[m_data] int start_q;
  array[len_q] int q_id;
  matrix[len_q,len_kappa_seq] q;
  matrix[len_q,len_kappa_seq] nu_main;
  matrix[len_q,len_kappa_seq] nu_urban;
  
  vector[m_data] y;                     // direct estimates
  vector[m_data] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  //constants for Satt approximation
  matrix[m_data,2] Cons;
  
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
  real<lower=-1> kappa; // urban contribution to variance of individuals
  vector<lower=0>[2] sig_u;
  //vector<lower=0>[2] phi;
  real<lower=0,upper=1> phi;
}

transformed parameters {
  vector[m] theta;
  vector[m] theta_drop;
  vector[m_data] theta_data;
  vector[m] log_sig2;
  vector[m_data] log_sig2_data;
  vector<lower=0>[m_data] v_raw;
  vector<lower=0>[m_data] df;
  vector[len_q] nu;
  vector[len_q] delta;
  real log_kappa_shift;
  
  vector[m] b1;
  vector[m] b2;
  
  log_kappa_shift = log(kappa+1);
  
  // random effects:
  b1 = sig_u[1]*(sqrt(phi/car_scale)*s1 + sqrt(1-phi)*u1);
  b2 = sig_u[2]*u2;
  
  //b1 = sig_u[1]*(sqrt(phi[1]/car_scale)*s1 + sqrt(1-phi[1])*u1);
  //b2 = sig_u[2]*(sqrt(phi[2]/car_scale)*s2 + sqrt(1-phi[2])*u2);

  // Mean and variance models
  theta = X * beta + b1;
  log_sig2 = Z * gamma + b2;
  
  // take out strata effect
  theta_drop = theta - X[,2]*beta[2];
  
  for(i in 1:m_data){
   vector[q_per_area[i]] q_tmp;
   vector[q_per_area[i]] nu_tmp;
   vector[q_per_area[i]] delta_tmp;
    
    theta_data[i] = theta[data_areas[i]];
    log_sig2_data[i] = log_sig2[data_areas[i]];
    
    q_tmp = interpolate_vec(kappa, kappa_seq, q[start_q[i]:(start_q[i] + q_per_area[i] -1),]);
    nu_tmp = theta_drop[data_areas[i]]*interpolate_vec(kappa, kappa_seq, nu_main[start_q[i]:(start_q[i] + q_per_area[i] -1),]) + beta[2]*interpolate_vec(kappa, kappa_seq, nu_urban[start_q[i]:(start_q[i] + q_per_area[i] -1),]);
    delta_tmp = square(nu_tmp./q_tmp);
    
    v_raw[i] = sum(q_tmp.*(1+delta_tmp))/sum(q_tmp.*q_tmp.*(1+2*delta_tmp))*v_hat_scaled[i]/exp(log_sig2_data[i]);
    df[i] = square(sum(q_tmp.*(1+delta_tmp)))/sum(q_tmp.*q_tmp.*(1+2*delta_tmp));
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
  
  log_kappa_shift ~ normal(0,0.5);
  
  beta ~ normal(0, 2);
  gamma[1] ~ normal(-5, 1);   // intercept-like prior; adjust if Z includes intercept
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,1);
  }
 
  v_raw ~ gamma(0.5*df,0.5);
  //y ~ normal(theta_data,sqrt(exp(log_sig2_data).*Cons));
  y ~ normal(theta_data,sqrt(exp(log_sig2_data).*(Cons[,1]+ (1+kappa)*Cons[,2])));
 
}

