data{
  int<lower=1> n_obs;
  int<lower=1> n_clr;
  int<lower=1, upper=n_clr> obs2clr[n_obs];
  vector[n_obs] y;
  // vector[n_obs] y_sd;
  vector[n_obs] x;
  vector[n_obs] x_sd;
}

parameters{
  // vector[n_obs] y_lat;
  vector[n_obs] x_lat;
  vector[n_clr] b0;
  vector[n_clr] b1;
  vector[n_clr] b2;
  vector<lower=0>[n_clr] sigma;
}

transformed parameters{
  vector[n_obs] mu;
  for (i in 1:n_obs) mu[i] = b0[obs2clr[i]] + b1[obs2clr[i]] * x_lat[i] + b2[obs2clr[i]] * (x_lat[i])^2;
}

model{
  x ~ normal(x_lat, x_sd);
  // y ~ normal(y_lat, y_sd);    
  for (i in 1:n_obs) y[i] ~ normal(mu[i], sigma[obs2clr[i]]);
  
  b0 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,10);
  // sigma ~ gamma(5,1); 
  sigma ~ normal(0,100);
  x_lat ~ normal(0, 50);
}
