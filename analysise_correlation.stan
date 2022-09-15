data {
	int<lower=1> n;
	int<lower=1> k;
	vector[k] x[n];
	vector[n] x_obs;
	vector[6] x_sd[n];	
}

parameters {
	vector[k] mu;
	vector[6] x_dash[n];
   vector<lower=0>[k] sigma;
	real<lower=1> nu; 
	cholesky_factor_corr[k] Lcorr;// cholesky factor (L_u matrix for R)
}

transformed parameters {
   
	corr_matrix[k] R; // correlation matrix
	cov_matrix[k] cov; // VCV matrix
	R = multiply_lower_tri_self_transpose(Lcorr);
	cov = quad_form_diag(R, sigma); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)
}

model {
   for(j in 1:6){
      for (i in 1:n){
         if (!x_obs[i] && j==1)
         	x_dash[i,j] ~ normal(0,100);
         else 
				x_dash[i,j] ~ normal(x[i,j], x_sd[i,j]);
      }
   }
	for (i in 1:n){
		if (x_obs[i]) // no missing data
			[x_dash[i,1], x_dash[i,2], x_dash[i,3],
         x_dash[i,4], x_dash[i,5], x_dash[i,6], 
         x[i,7], x[i,8], x[i,9]] ~ multi_student_t(nu, mu, cov);
		else 
         [x_dash[i,2], x_dash[i,3],
         x_dash[i,4], x_dash[i,5], x_dash[i,6], 
         x[i,7], x[i,8], x[i,9]] ~ multi_student_t(nu, mu[2:k], cov[2:k,2:k]);
	}
	
	// x ~ multi_normal(mu, cov);
	// sigma ~ normal(0,1);
	Lcorr ~ lkj_corr_cholesky(2);
	// mu ~ normal(0,1);
	nu ~ gamma(2,.1);

}
