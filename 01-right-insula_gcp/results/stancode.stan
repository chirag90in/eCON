// generated with brms 2.12.0
functions {

  /* compute the logm1 link 
   * Args: 
   *   p: a positive scalar
   * Returns: 
   *   a scalar in (-Inf, Inf)
   */ 
   real logm1(real y) { 
     return log(y - 1);
   }
  /* compute the inverse of the logm1 link 
   * Args: 
   *   y: a scalar in (-Inf, Inf)
   * Returns: 
   *   a positive scalar
   */ 
   real expp1(real y) { 
     return exp(y) + 1;
   }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  vector[N] Z_1_3;
  vector[N] Z_1_4;
  vector[N] Z_1_5;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  vector[N] Z_3_3;
  vector[N] Z_3_4;
  vector[N] Z_3_5;
  int<lower=1> NC_3;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // residual SD
  real<lower=1> nu;  // degrees of freedom or shape
  // parameters for student-t distributed group-level effects
  real<lower=1> df_1;
  vector<lower=0>[N_1] udf_1;
  real<lower=1> df_2;
  vector<lower=0>[N_2] udf_2;
  real<lower=1> df_3;
  vector<lower=0>[N_3] udf_3;
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // standardized group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
}
transformed parameters {
  vector[N_1] dfm_1;
  vector[N_2] dfm_2;
  vector[N_3] dfm_3;
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  vector[N_1] r_1_3;
  vector[N_1] r_1_4;
  vector[N_1] r_1_5;
  vector[N_2] r_2_1;  // actual group-level effects
  matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_1;
  vector[N_3] r_3_2;
  vector[N_3] r_3_3;
  vector[N_3] r_3_4;
  vector[N_3] r_3_5;
  dfm_1 = sqrt(df_1 * udf_1);
  dfm_2 = sqrt(df_2 * udf_2);
  dfm_3 = sqrt(df_3 * udf_3);
  // compute actual group-level effects
  r_1 = rep_matrix(dfm_1, M_1) .* (diag_pre_multiply(sd_1, L_1) * z_1)';
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  r_1_3 = r_1[, 3];
  r_1_4 = r_1[, 4];
  r_1_5 = r_1[, 5];
  r_2_1 = dfm_2 .* (sd_2[1] * (z_2[1]));
  // compute actual group-level effects
  r_3 = rep_matrix(dfm_3, M_3) .* (diag_pre_multiply(sd_3, L_3) * z_3)';
  r_3_1 = r_3[, 1];
  r_3_2 = r_3[, 2];
  r_3_3 = r_3[, 3];
  r_3_4 = r_3[, 4];
  r_3_5 = r_3[, 5];
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n] + r_1_4[J_1[n]] * Z_1_4[n] + r_1_5[J_1[n]] * Z_1_5[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_3_2[J_3[n]] * Z_3_2[n] + r_3_3[J_3[n]] * Z_3_3[n] + r_3_4[J_3[n]] * Z_3_4[n] + r_3_5[J_3[n]] * Z_3_5[n];
  }
  // priors including all constants
  target += student_t_lpdf(b | 3, 0, 10);
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += gamma_lpdf(nu | 3.325,0.1)
    - 1 * gamma_lccdf(1 | 3.325,0.1);
  target += gamma_lpdf(df_1 | 3.325,0.1);
  target += inv_chi_square_lpdf(udf_1 | df_1);
  target += gamma_lpdf(df_2 | 3.325,0.1);
  target += inv_chi_square_lpdf(udf_2 | df_2);
  target += gamma_lpdf(df_3 | 3.325,0.1);
  target += inv_chi_square_lpdf(udf_3 | df_3);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 5 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 2);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_2[1] | 0, 1);
  target += student_t_lpdf(sd_3 | 3, 0, 10)
    - 5 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(to_vector(z_3) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_3 | 2);
  // likelihood including all constants
  if (!prior_only) {
    target += student_t_lpdf(Y | nu, mu, sigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
}

