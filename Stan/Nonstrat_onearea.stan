
data {
  real y;                     // direct estimates
  real v_hat_scaled;         // variance estimates, rescaled for chi square approximation
  
  //constants 
  real df;
  real Cons;
}

// The parameters accepted by the model. 
parameters {
  real theta;
  real log_sig2;
}

transformed parameters {
  real<lower=0> v;
  real<lower=0> v_raw;

  v = exp(log_sig2)*Cons;
  v_raw = v_hat_scaled/v;
 
}

// The model to be estimated.
model {
  
  theta ~ normal(0, 2);
  log_sig2 ~ normal(0.5, 0.5);
  
  v_raw ~ gamma(0.5*df,0.5);
  y ~ normal(theta,sqrt(v));

}
