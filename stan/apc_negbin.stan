// Modele age-periode-cohorte sous loi binomiale negative.
// Identification : kappa[1] = 0, cohortes extremes fixees a zero,
// et effet de cohorte orthogonal au niveau et a la tendance lineaire.

data {
  int<lower=2> A;
  int<lower=2> T;
  int<lower=4> C;
  array[A, T] int<lower=0> D;
  matrix[A, T] log_E;
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
  vector[A] alpha;
  real drift;
  real<lower=0> sigma_kappa;
  vector[T - 1] kappa_free;

  vector[C - 4] gamma_coord;
  real<lower=-0.99, upper=0.99> pacf1;
  real<lower=-0.99, upper=0.99> pacf2;
  real<lower=0> sigma_gamma;
  real<lower=0> inv_phi;
}

transformed parameters {
  vector[T] kappa;
  vector[C] gamma;
  real psi1 = pacf1 * (1 - pacf2);
  real psi2 = pacf2;
  real<lower=0> phi = inv(inv_phi);

  kappa[1] = 0;
  kappa[2:T] = kappa_free;

  gamma[1] = 0;
  gamma[2:(C - 1)] = Q_cohort * gamma_coord;
  gamma[C] = 0;
}

model {
  vector[A * T] eta;

  alpha ~ normal(-5, 3);
  drift ~ normal(0, 0.1);
  sigma_kappa ~ exponential(10);
  for (t in 2:T) {
    kappa[t] ~ normal(kappa[t - 1] + drift, sigma_kappa);
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
      eta[i] = alpha[a] + kappa[t] + gamma[cohort_idx[a, t]];
    }
  }
  D_long ~ neg_binomial_2_log(log_E_long + eta, phi);
}
