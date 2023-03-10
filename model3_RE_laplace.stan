//
// This Stan program defines a simple random effect model with multivariate outcome, 
// 
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; //number of individual
  int<lower=0> n_complete; // number of complete cases
  int<lower=0> n_missing; // number of missing cells - among the incompleted cases
  int<lower=0> n_obs; // number of observed cells - among the incompleted case
  int<lower =0> K;
  vector[K] y_all[N];
  matrix[K,K] x_mat[N];
  
  // infor regarding missingness
  int missingRow[n_missing]; // row indices of missing items in the full data
  int missingCol[n_missing]; // column indices of missing items in the full data
  int obsCol[n_obs] ; // column indices of observed items for partially obsereved observations
  int obsRow[n_obs] ; // row indices of observed items for partially obsereved observations
  
  // prior setup
  real sigma_prior_intercept;
  real sigma_prior_logsigma;
  real sigma_prior_logsigmaAlpha;
  real sigma_prior_gamma;
  
  real pscore[N]; // propensity score
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[K] beta; // per-outcome effect
  vector[K] gamma; 
  vector[K] intercept; // intercept
  vector[N] alpha_norm; // individual random effect
  real logsigma_alpha;
  real logsigma[K];// variance, for now assume the outcome are independent
  //real <lower = 0> sigma;
  real ymiss[n_missing];
}
transformed parameters{
  vector[N] alpha;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma[K]; // variance, for now assume the outcome are
  
  sigma_alpha = exp(logsigma_alpha);
  alpha = sigma_alpha*alpha_norm;
  sigma = exp(logsigma);
}

model {
  vector[K] fullData[N]; // complete data
  //matrix[K, K] L_Sigma;
  vector[K] mean_y[N];
  for(i in 1:N){
    mean_y[i] = intercept + alpha[i] + x_mat[i]*beta+ pscore[i]*gamma;
  }
  
  
  intercept ~ normal(0,sigma_prior_intercept);
 
  //sigma_alpha ~ cauchy(0, 5);
  //sigma ~ cauchy(0,5);
  logsigma_alpha~ normal(0,sigma_prior_logsigmaAlpha);
  logsigma~ normal(0,sigma_prior_logsigma);
  gamma ~ normal(0,sigma_prior_gamma);
  //sigma_total ~ cauchy(0,5);

  beta ~ double_exponential(0, 1);
  
  alpha_norm ~ normal(0, 1);
  
  fullData[1:n_complete] = y_all[1:n_complete];
  if(n_missing>0){
    for(r in 1:n_missing){
    fullData[missingRow[r],missingCol[r]]  = ymiss[r];
  }
  for(r in 1:n_obs){
    fullData[obsRow[r],obsCol[r]]  = y_all[obsRow[r],obsCol[r]];
  }
  }
  
  
  for (i in 1:N){
    fullData[i] ~ normal( mean_y[i],sigma);
    //y_all[i] ~ normal( mean_y[i],sigma);
  }
}
