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
  vector[m_data] v_hat;         // variance estimates, rescaled for chi square approximation
  
  // Covariates for the mean model
  int<lower=1> p_mean;             // number of mean covariates
  matrix[m, p_mean] X;             // covariate matrix for all areas

  int<lower=0,upper=1> bym2_mean;   // 1 = BYM2 for mean model, 0 = IID
  
   // Adjacency (undirected) for ICAR admin2
  int<lower=0> N_edges;            // number of edges
  array[N_edges] int<lower=1, upper=m> node_1;
  array[N_edges] int<lower=1, upper=m> node_2;
  real<lower=0> car_scale;
}

// The parameters accepted by the model. 
parameters {
  vector[p_mean] beta;             // coefficients for mean model
  vector[m] u1;
  vector[m] z_s1;   // unconstrained ICAR
  real<lower=0> sig_u;
  real<lower=0,upper=1> phi1;
}

transformed parameters {
  vector[m] theta;
  vector[m_data] theta_data;
  vector[m] s1; // ICAR
  vector[m] s1_raw; // ICAR
  vector[m] b1;
  
  s1_raw = z_s1 - mean(z_s1);      // sum-to-zero enforced here

  // random effects:
  if(bym2_mean==1){
    s1 = s1_raw / sqrt(car_scale);   // BYM2 standardization
    b1 = sig_u*(sqrt(phi1)*s1 + sqrt(1-phi1)*u1);
  }else{
    b1 = sig_u*u1;
  }
  
  theta = X * beta + b1;
  
  for(a in 1:m_data){
    theta_data[a] = theta[data_areas[a]];
  }
 
}

// The model to be estimated.
model {
  target += pc_sigma_lp(sig_u, 1, 0.01);
  phi1 ~ beta(0.5,1);

  u1 ~ normal(0,1);
  z_s1 ~ normal(0, 1);

  if(bym2_mean==1){
    target += icar_lp(s1_raw, node_1, node_2);
  }
  
  beta ~ normal(0, 2);
  y ~ normal(theta_data,sqrt(v_hat));

}
