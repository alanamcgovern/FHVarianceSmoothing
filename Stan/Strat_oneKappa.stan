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
  int m;
  int m_data;
  array[m_data] int data_areas;
  vector[m_data] y;                     // direct estimates
  vector[m_data] v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  //constants for Satt approximation
  int lenq;
  int len_kappa_seq;
  vector[len_kappa_seq] kappa_seq;
  array[lenq] int q_id;
  matrix[m_data,2] Cons;
  array[m_data] int q_start;
  array[m_data] int q_per_area;
  matrix[lenq,len_kappa_seq] q_mat;
  matrix[lenq,len_kappa_seq] nu;
 //vector[lenq] q_mat;
 //vector[lenq] nu;
  
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
  vector[m] theta;
  vector[m] log_v;
  vector[m_data] theta_data;
  vector<lower=0>[m_data] v_data;
  vector<lower=0>[m_data] v_raw;
  vector[m] s1; // ICAR
  vector[m] s2;
  vector[m] b1;
  vector[m] b2;
  vector[lenq] q;
  vector[lenq] delta;
  matrix[m_data,2] sums;
  vector[m_data] df;
  real kappa_trans;
  
  kappa_trans = log(kappa+1);
  
  s1 = (s1_raw - mean(s1_raw))/sqrt(car_scale);
  s2 = (s2_raw - mean(s2_raw))/sqrt(car_scale);
  
  
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
  log_v = Z * gamma + b2;
  
  //q = q_mat;
  //delta = square(nu*beta[2])./exp(log_v[q_id]);
  q = interpolate_vec(kappa,kappa_seq,q_mat);
  delta = square(interpolate_vec(kappa,kappa_seq,nu)*beta[2])./exp(log_v[q_id]);

  for(a in 1:m_data){
    theta_data[a] = theta[data_areas[a]];
    v_data[a] = exp(log_v[data_areas[a]]);

    sums[a,1] = sum(q[q_start[a]:(q_start[a]+q_per_area[a]-1)].*(1+delta[q_start[a]:(q_start[a]+q_per_area[a]-1)]));
    sums[a,2] = sum(square(q[q_start[a]:(q_start[a]+q_per_area[a]-1)]).*(1+2*delta[q_start[a]:(q_start[a]+q_per_area[a]-1)]));
  }
  
  df = square(sums[,1])./sums[,2];
  v_raw = v_hat_scaled./v_data.*sums[,1]./sums[,2].*(Cons[,1]+kappa*Cons[,2]);

}

// The model to be estimated.
model {
  target += pc_sigma_lp(sig_u[1], 1, 0.01);
  target += pc_sigma_lp(sig_u[2], 1, 0.01);
  
  phi1 ~ beta(0.25,1);
  phi2 ~ beta(0.25,1);
  
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
  
  gamma[1] ~ normal(-3, 0.5);
  if(p_var>1){
     gamma[2:p_var] ~ normal(0,1);
  }
 
  v_raw ~ gamma(0.5*df,0.5);
  y ~ normal(theta_data,sqrt(v_data));

}
