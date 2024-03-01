data {
  real<lower=0> beta_a; 
  real<lower=0> beta_b; 
  int<lower=0> n_total_1; 
  int<lower=0> n_total_2;
  int<lower=0> n_hits_1;
  int<lower=0> n_hits_2; 
  int<lower=0, upper=1> sample_prior_only;
}
parameters { 
  real<lower=0, upper=1> theta_1;
  real<lower=0, upper=1> theta_2;
}
model {
  theta_1 ~ beta(beta_a, beta_b);
  theta_2 ~ beta(beta_a, beta_b);
  
  if (sample_prior_only != 1) {
    n_hits_1 ~ binomial(n_total_1, theta_1);
    n_hits_2 ~ binomial(n_total_2, theta_2);
  }
}
generated quantities {
  real diff;
  real OR;
  diff = theta_1 - theta_2;
  OR = (theta_1 / (1 - theta_1)) / (theta_2 / (1 - theta_2));
}
