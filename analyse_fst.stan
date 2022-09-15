data {
  int<lower=1> n_obs;
  int<lower=1> n_clr;
  int<lower=1> n_eye;
  int<lower=1> n_pat;
  int<lower=1, upper=n_eye> obs2eye[n_obs];
  int<lower=1, upper=n_pat> obs2pat[n_obs];
  int<lower=1, upper=n_clr> obs2col[n_obs];
  int<lower=1, upper=n_pat> eye2pat[n_eye];
  real x[n_obs];
  int<lower=0, upper=1> y[n_obs];
}

parameters {
  // real a0;
  real beta;
  real a_eye[n_eye];
  real<lower=0> sigma_eye; 
  real a_pat[n_pat];
  real<lower=0> sigma_pat;
  real a_clr[n_clr];
  real a_eyeXclr[n_eye, n_clr];
  real<lower=0> simga_eyeXclr;
  real<lower=0, upper=1> guess[n_eye];
}

model {
  for (i in 1:n_obs){
     y[i] ~ bernoulli(
        .5*guess[obs2eye[i]] + (1-guess[obs2eye[i]]) * inv_logit(
          a_clr[obs2col[i]] +
          a_pat[obs2pat[i]] +
          a_eye[obs2eye[i]] +
          a_eyeXclr[obs2eye[i], obs2col[i]] +
          beta * x[i] )
        );
  }
  a_pat ~ normal(0, sigma_pat);
  sigma_pat ~ gamma(1.64,.32);
  a_eye ~ normal(0, sigma_eye);
  sigma_eye ~ gamma(1.64, .32);
  a_clr ~ student_t(3,0,1);
    // interactions
  for(i1 in 1:n_eye){ //index line
    for(i2 in 1:n_clr){ //index cond and topo
      a_eyeXclr[i1,i2] ~ normal(0, simga_eyeXclr);
    }
  }
  simga_eyeXclr ~ gamma(1.64, .32);
  guess ~ beta(1,9);
}


generated quantities {
  real b0=0.0;
  real b_eye[n_eye];
  real b_pat[n_pat];
  real b_clr[n_clr];
  real b_eyeXclr[n_eye, n_clr];
  real y_pred[n_eye, n_clr];
  
  {
    real m[n_eye, n_pat, n_clr];
    for(e in 1:n_eye){
      for(p in 1:n_pat)
        for(c in 1:n_clr){
          m[e,p,c] = a_eye[e] + a_pat[p] + a_clr[c] + a_eyeXclr[e,c];
          b0 += m[e,p,c];
        }
    }
  
    for (e in 1:n_eye){
      for (c in 1:n_clr)
      y_pred[e,c] = a_eye[e] + a_pat[eye2pat[e]] + a_clr[c] + a_eyeXclr[e,c];
    }
    b0 = b0 / num_elements(m);
    for(e in 1:n_eye) b_eye[e] = mean(to_array_1d(m[e,:,:])) - b0; 
    for(p in 1:n_pat) b_pat[p] = mean(to_array_1d(m[:,p,:])) - b0;   
    for(c in 1:n_clr) b_clr[c] = mean(to_array_1d(m[:,:,c])) - b0; 
    for (e in 1:n_eye){
      for (c in 1:n_clr){
        b_eyeXclr[e,c] = mean(to_array_1d(m[e,:,c])) - (b_eye[e] + b_clr[c] + b0);
      }
    }
  }
}


