data {
  int<lower=1> n1;
  int<lower=1> n2;
  int<lower=1> n_eye;
  int<lower=1, upper=n_eye> obs2eye1[n1];
  int<lower=1, upper=n_eye> obs2eye2[n2];
  int<lower=1> n_pat;
  int<lower=1, upper=n_pat> obs2pat1[n1];
  int<lower=1, upper=n_pat> obs2pat2[n2];
  int<lower=1, upper=n_pat> eye2pat[n_eye];
  int<lower=1, upper=n_eye> eye2va[n_eye];
  matrix<lower=0, upper=1>[n_eye,2] x_obs;
  real x1[n1];
  real x2_mu[n2];
  real x2_sd[n2];
}

parameters {
  vector[n_eye] mu_eye_x1;
  vector[n_eye] mu_eye_x2;
  vector[n_pat] mu_pat_x1;
  vector[n_pat] mu_pat_x2;
  // real mu_pat_x1_mu0;
  // real<lower=0> mu_pat_x1_sd;
  // real mu_eye_x1_mu0;
  // real<lower=0> mu_eye_x1_sd;
  // real mu_pat_x2_mu0;
  // real<lower=0> mu_pat_x2_sd;
  // real mu_eye_x2_mu0;
  // real<lower=0> mu_eye_x2_sd;  
  real<lower=0> sd_x1;
  vector[2] mu;
  vector<lower=0>[2] sigma;
  real<lower=1> nu; 
  cholesky_factor_corr[2] Lcorr[4];// cholesky factor (L_u matrix for R)
}

transformed parameters {
  vector[n_eye] mu_x1;
  vector[n_eye] mu_x2;
  corr_matrix[2] R[4]; // correlation matrix
  cov_matrix[2] cov[4]; // VCV matrix
  
  for (eye in 1:n_eye){
    mu_x1[eye] = mu_eye_x1[eye] + mu_pat_x1[eye2pat[eye]];
    mu_x2[eye] = mu_eye_x2[eye] + mu_pat_x2[eye2pat[eye]];
  }

  for (v in 1:4){
    R[v] = multiply_lower_tri_self_transpose(Lcorr[v]);
    cov[v] = quad_form_diag(R[v], sigma); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)    
  }

}

model {
  for (i in 1:n1){
    x1[i] ~ normal(mu_x1[obs2eye1[i]], sd_x1);
  }

  for (i in 1:n2){
    x2_mu[i] ~ normal(mu_x2[obs2eye2[i]], x2_sd[i]);
  }

  for (eye in 1:n_eye){
    if (x_obs[eye, 1] && x_obs[eye, 2])
      [mu_x1[eye], mu_x2[eye]] ~ multi_student_t(nu, mu, cov[eye2va[eye]]);
    else if (x_obs[eye, 1])
      mu_x1[eye] ~ student_t(nu, mu[1], cov[eye2va[eye]][1,1]);
    else if (x_obs[eye, 2])
      mu_x2[eye] ~ student_t(nu, mu[2], cov[eye2va[eye]][2,2]);
  }

  
  // retina thickness
  sd_x1 ~ normal(0,10);
  mu_pat_x1 ~ normal(200,50);
  mu_eye_x1 ~ normal(0,20);  
  
  // FST
  mu_pat_x2 ~ normal(0,30);
  mu_eye_x2 ~ normal(0,10);
 
  sigma[1] ~ normal(0,200);
  sigma[2] ~ normal(0,20);
  for (v in 1:4){
    Lcorr[v] ~ lkj_corr_cholesky(2);    
  }

  nu ~ gamma(2,.1);
}

generated quantities{
  vector[2] pred[4];

  for (v in 1:4){
    pred[v] = multi_student_t_rng(nu, mu, cov[v]);
  }
}
