data {
  int<lower=1> n_obs;
  // int<lower=2> k; n_var = 3 before, peak, after
  matrix[n_obs,3] x;
  matrix<lower=0, upper=1>[n_obs,3] x_obs;
  int<lower=1> n_stp;
  int<lower=1, upper=n_stp> obs2stp[n_obs];
  int<lower=1> n_eye;
  int<lower=1, upper=n_eye> obs2eye[n_obs];
  int<lower=1> n_pat;
  int<lower=1, upper=n_pat> obs2pat[n_obs];
  int<lower=1, upper=n_pat> eye2pat[n_eye];
}

parameters {
  real m_eye[n_eye,3];
  real m_pat[n_pat,3];
  real m_stp3[3,3]; // for stp 1,2,3 model before, peak, after
  real m_stp2[2,2]; // for stp 4,5 model before, after
  real m_eyeXstp3[n_eye,3,3]; // for stp 1,2,3 model before, peak, after
  real m_eyeXstp2[n_eye,2,2]; // for stp 4,5 model before, after
  real<lower=0> s_eye; 
  real<lower=0> s_pat;
  real<lower=0> s_stp;
  real<lower=0> s_eyeXstp; 
  vector<lower=0>[3] sigma;
  real<lower=1> nu; 
  cholesky_factor_corr[3] Lcorr3[3];// cholesky factor (L_u matrix for R)
  cholesky_factor_corr[2] Lcorr2[2];
}

transformed parameters {
  vector[3] mu[n_obs];
  corr_matrix[3] R3[3]; // correlation matrix for stp 1,2,3
  corr_matrix[2] R2[2]; // correlation matrix for stp 4,5
  cov_matrix[3] cov3[3]; // VCV matrix
  cov_matrix[2] cov2[2];
  cov_matrix[2] cov12[3];
  cov_matrix[2] cov13[3];
  cov_matrix[2] cov23[3];
  
  for (stp in 1:5){
    if (stp<4){
      R3[stp] = multiply_lower_tri_self_transpose(Lcorr3[stp]);
      cov3[stp] = quad_form_diag(R3[stp], sigma); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)
      // cov12
      cov12[stp][1,1] = cov3[stp][1,1];
      cov12[stp][1,2] = cov3[stp][1,2];
      cov12[stp][2,1] = cov3[stp][2,1];
      cov12[stp][2,2] = cov3[stp][2,2];
      // cov13
      cov13[stp][1,1] = cov3[stp][1,1];
      cov13[stp][1,2] = cov3[stp][1,3];
      cov13[stp][2,1] = cov3[stp][3,1];
      cov13[stp][2,2] = cov3[stp][3,3];
      // cov23
      cov23[stp][1,1] = cov3[stp][2,2];
      cov23[stp][1,2] = cov3[stp][2,3];
      cov23[stp][2,1] = cov3[stp][3,2];
      cov23[stp][2,2] = cov3[stp][3,3];      
    } else if (stp>3){
      R2[stp-3] = multiply_lower_tri_self_transpose(Lcorr2[stp-3]);
      cov2[stp-3] = quad_form_diag(R2[stp-3], [sigma[1], sigma[3]]); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)
    }
  }

  for(i in 1:n_obs){
    for (k in 1:3) mu[i,k] = m_eye[obs2eye[i],k] + m_pat[obs2pat[i],k];
    if (obs2stp[i]<4){
      for (k in 1:3){
        mu[i,k] = mu[i,k] +
          m_stp3[obs2stp[i],k] + 
          m_eyeXstp3[obs2eye[i], obs2stp[i],k];               
      }
    } else if (obs2stp[i]>3){
      mu[i,1] = mu[i,1] +
        m_stp2[obs2stp[i]-3,1] + 
        m_eyeXstp2[obs2eye[i], obs2stp[i]-3,1];
      mu[i,3] = mu[i,3] +
        m_stp2[obs2stp[i]-3,2] + 
        m_eyeXstp2[obs2eye[i], obs2stp[i]-3,2];  
      mu[i,2] = 0.0; // this is undefined, do I need to include this?
    }
  }
}

model {
  for (i in 1:n_obs){
    if (obs2stp[i]<4){
      if (x_obs[i,1] && x_obs[i,2] && x_obs[i,3]){
        x[i] ~ multi_student_t(nu, mu[i], cov3[obs2stp[i]]);      
      } else if (x_obs[i,1] && x_obs[i,2]){
        [x[i,1], x[i,2]] ~ multi_student_t(nu, [mu[i,1],mu[i,2]], cov12[obs2stp[i]]);
      } else if (x_obs[i,2] && x_obs[i,3]){
        [x[i,2], x[i,3]]  ~ multi_student_t(nu, [mu[i,2],mu[i,3]], cov23[obs2stp[i]]);
      } else if (x_obs[i,1] && x_obs[i,3]){
        [x[i,1], x[i,3]]  ~ multi_student_t(nu, [mu[i,1],mu[i,3]], cov13[obs2stp[i]]);      
      } else if (x_obs[i,1]){
        x[i,1] ~ student_t(nu, mu[i,1], sqrt(cov3[obs2stp[i]][1,1]));
      } else if (x_obs[i,2]){
        x[i,2] ~ student_t(nu, mu[i,2], sqrt(cov3[obs2stp[i]][2,2]));
      } else if (x_obs[i,3]){
        x[i,3] ~ student_t(nu, mu[i,3], sqrt(cov3[obs2stp[i]][3,3]));
      }      
    } else if (obs2stp[i]>3){
      if (x_obs[i,1] && x_obs[i,3]){
        [x[i,1], x[i,3]] ~ multi_student_t(nu, [mu[i,1],mu[i,3]], cov2[obs2stp[i]-3]);
      } else if (x_obs[i,1]){
        x[i,1] ~ student_t(nu, mu[i,1], sqrt(cov2[obs2stp[i]-3][1,1]));
      } else if (x_obs[i,3]){
        x[i,3] ~ student_t(nu, mu[i,3], sqrt(cov2[obs2stp[i]-3][2,2]));
      }
    }
  }

  sigma ~ normal(0,5);
  for (stp in 1:n_stp){
    if (stp<4){
      Lcorr3[stp] ~ lkj_corr_cholesky(2);
      m_stp3[stp,] ~ normal(0,s_stp);
      for (eye in 1:n_eye) m_eyeXstp3[eye,stp,] ~ normal(0, s_eyeXstp);
    } else if (stp>3){
      Lcorr2[stp-3] ~ lkj_corr_cholesky(2);
      m_stp2[stp-3,] ~ normal(0,s_stp);
      for (eye in 1:n_eye) m_eyeXstp2[eye,stp-3,] ~ normal(0, s_eyeXstp);
    }
  }
  for (eye in 1:n_eye) m_eye[eye,] ~ normal(0,s_eye);
  for (pat in 1:n_pat) m_pat[pat,] ~ normal(0,s_pat);    
  s_pat ~ normal(0,10);
  s_stp ~ normal(0,10);
  s_eye ~ normal(0,10);
  s_eyeXstp ~ normal(0,10);
  nu ~ gamma(2,.1);

}


generated quantities {
  real b0[3];
  real b_eye[n_eye,3];
  real b_pat[n_pat,3];
  real b_stp3[3,3];
  real b_stp2[2,2];
  real b_eyeXstp3[n_eye, 3, 3];
  real b_eyeXstp2[n_eye, 2, 2];
  real pred3[n_eye, 3, 3];
  real pred2[n_eye, 2, 2];
  real m3[n_eye, n_pat, 3, 3];
  real m2[n_eye, n_pat, 2, 2];

  for (k in 1:3) b0[k] = 0.0;
  for(eye in 1:n_eye){
    for(stp in 1:n_stp){
      if (stp<4){
        for (k in 1:3){
          pred3[eye,stp,k] = m_eye[eye,k] + m_pat[eye2pat[eye],k] + m_stp3[stp,k] + m_eyeXstp3[eye,stp,k];  
        }
      } else if (stp>3){
        pred2[eye,stp-3,1] = m_eye[eye,1] + m_pat[eye2pat[eye],1] + m_stp2[stp-3,1] + m_eyeXstp2[eye,stp-3,1];  
        pred2[eye,stp-3,2] = m_eye[eye,3] + m_pat[eye2pat[eye],3] + m_stp2[stp-3,2] + m_eyeXstp2[eye,stp-3,2];  
      }
      for(pat in 1:n_pat){
        if (stp<4){
          for (k in 1:3){
            m3[eye,pat,stp,k] = m_eye[eye,k] + m_pat[pat,k] + m_stp3[stp,k] + m_eyeXstp3[eye,stp,k];
            b0[k] += m3[eye,pat,stp,k];
          }
        }
        if (stp>3){
          m2[eye,pat,stp-3,1] = m_eye[eye,1] + m_pat[pat,1] + m_stp2[stp-3,1] + m_eyeXstp2[eye,stp-3,1];
          m2[eye,pat,stp-3,2] = m_eye[eye,3] + m_pat[pat,3] + m_stp2[stp-3,2] + m_eyeXstp2[eye,stp-3,2];
          b0[1] += m2[eye,pat,stp-3,1];
          b0[3] += m2[eye,pat,stp-3,2];
        }                    
      }
    }
  }

    
  b0[1] = b0[1] / (num_elements(m3[:,:,:,1]) + num_elements(m2[:,:,:,1]));
  b0[2] = b0[2] / num_elements(m3[:,:,:,2]);  
  b0[3] = b0[3] / (num_elements(m3[:,:,:,3]) + num_elements(m2[:,:,:,2]));

  for (eye in 1:n_eye){
    b_eye[eye,1] = mean(append_array(to_array_1d(m3[eye,:,:,1]), to_array_1d(m2[eye,:,:,1]))) - b0[1];
    b_eye[eye,2] = mean(to_array_1d(m3[eye,:,:,2])) - b0[2];
    b_eye[eye,3] = mean(append_array(to_array_1d(m3[eye,:,:,3]), to_array_1d(m2[eye,:,:,2]))) - b0[3];
  }

  for(pat in 1:n_pat) {
    b_pat[pat,1] = mean(append_array(to_array_1d(m3[:,pat,:,1]), to_array_1d(m2[:,pat,:,1]))) - b0[1];
    b_pat[pat,2] = mean(to_array_1d(m3[:,pat,:,2])) - b0[2];
    b_pat[pat,3] = mean(append_array(to_array_1d(m3[:,pat,:,3]), to_array_1d(m2[:,pat,:,2]))) - b0[3];
  }

  for(stp in 1:n_stp){
    if(stp<4){
      for (k in 1:3){
        b_stp3[stp,k] = mean(to_array_1d(m3[:,:,stp,k])) - b0[k];  
      }
    } else if (stp>3){
        b_stp2[stp-3,1] = mean(to_array_1d(m2[:,:,stp-3,1])) - b0[1];  
        b_stp2[stp-3,2] = mean(to_array_1d(m2[:,:,stp-3,2])) - b0[3];  
    }
  }


    
  for (eye in 1:n_eye){
    for (stp in 1:n_stp){
      if (stp<4) {
        for (k in 1:3){
          b_eyeXstp3[eye,stp,k] = mean(to_array_1d(m3[eye,:,stp,k])) - (b_eye[eye,k] + b_stp3[stp,k] + b0[k]);          
        }
      } else if (stp>3){
          b_eyeXstp2[eye,stp-3,1] = mean(to_array_1d(m2[eye,:,stp-3,1])) - (b_eye[eye,1] + b_stp2[stp-3,1] + b0[1]);
          b_eyeXstp2[eye,stp-3,2] = mean(to_array_1d(m2[eye,:,stp-3,2])) - (b_eye[eye,3] + b_stp2[stp-3,2] + b0[3]);
      }
    }
  }
}
