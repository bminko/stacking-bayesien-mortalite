// Modele M6 : structure CBD augmentee d'un effet de cohorte AR(2).
// L'effet cohorte est orthogonal au niveau et a la tendance lineaire.
// La pente CBD est echantillonnee par ecart-type d'age, puis restituee
// comme coefficient par annee d'age dans kappa[2, ].

data {
  int<lower=2> A;
  int<lower=2> T;
  int<lower=4> C;
  array[A, T] int<lower=0> D;
  matrix[A, T] log_E;
  vector[A] age_centered;
  real<lower=0> age_scale;
  array[A, T] int<lower=1, upper=C> cohort_idx;
}

transformed data {
  int N = A * T;
  array[N] int D_long;
  vector[N] log_E_long;
  matrix[C - 2, C - 2] basis_seed = rep_matrix(0.0, C - 2, C - 2);
  matrix[C - 2, C - 2] Q_full;
  matrix[C - 2, C - 4] Q_cohort;
  for (t in 1:T) {
    for (a in 1:A) {
      int i = (t - 1) * A + a;
      D_long[i] = D[a, t];
      log_E_long[i] = log_E[a, t];
    }
  }
  basis_seed[, 1] = rep_vector(1.0, C - 2);
  for (i in 1:(C - 2)) {
    basis_seed[i, 2] = i - 0.5 * (C - 1);
  }
  for (j in 3:(C - 2)) {
    basis_seed[j - 2, j] = 1;
  }
  Q_full = qr_thin_Q(basis_seed);
  Q_cohort = block(Q_full, 1, 3, C - 2, C - 4);
}

parameters {
  matrix[2, T] kappa_scaled;
  vector[2] drift_scaled;
  vector<lower=0>[2] sigma_kappa_scaled;
  cholesky_factor_corr[2] L_Omega;

  vector[C - 4] gamma_coord;
  real<lower=-0.99, upper=0.99> pacf1;
  real<lower=-0.99, upper=0.99> pacf2;
  real<lower=0> sigma_gamma;
  real<lower=0> inv_phi;
}

transformed parameters {
  matrix[2, T] kappa;
  vector[2] drift;
  vector<lower=0>[2] sigma_kappa;
  matrix[2, 2] L_Sigma =
    diag_pre_multiply(sigma_kappa_scaled, L_Omega);
  corr_matrix[2] Omega = multiply_lower_tri_self_transpose(L_Omega);
  vector[C] gamma;
  real psi1 = pacf1 * (1 - pacf2);
  real psi2 = pacf2;
  real<lower=0> phi = inv(inv_phi);

  drift[1] = drift_scaled[1];
  drift[2] = drift_scaled[2] / age_scale;
  sigma_kappa[1] = sigma_kappa_scaled[1];
  sigma_kappa[2] = sigma_kappa_scaled[2] / age_scale;

  for (t in 1:T) {
    kappa[1, t] = kappa_scaled[1, t];
    kappa[2, t] = kappa_scaled[2, t] / age_scale;
  }

  gamma[1] = 0;
  gamma[2:(C - 1)] = Q_cohort * gamma_coord;
  gamma[C] = 0;
}

model {
  vector[A * T] eta;

  kappa_scaled[1, 1] ~ normal(-4, 2);
  kappa_scaled[2, 1] ~ normal(0, 2);
  drift_scaled ~ normal(0, 0.1);
  sigma_kappa_scaled ~ exponential(10);
  L_Omega ~ lkj_corr_cholesky(2);
  for (t in 2:T) {
    kappa_scaled[, t] ~ multi_normal_cholesky(
      kappa_scaled[, t - 1] + drift_scaled,
      L_Sigma
    );
  }

  pacf1 ~ normal(0, 0.5);
  pacf2 ~ normal(0, 0.5);
  sigma_gamma ~ exponential(10);
  gamma[2] ~ normal(0, 1);
  for (c in 3:(C - 1)) {
    gamma[c] ~ normal(psi1 * gamma[c - 1]
                      + psi2 * gamma[c - 2], sigma_gamma);
  }

  inv_phi ~ normal(0, 1);

  for (t in 1:T) {
    for (a in 1:A) {
      int i = (t - 1) * A + a;
      eta[i] = kappa[1, t] + age_centered[a] * kappa[2, t]
               + gamma[cohort_idx[a, t]];
    }
  }
  D_long ~ neg_binomial_2_log(log_E_long + eta, phi);
}
