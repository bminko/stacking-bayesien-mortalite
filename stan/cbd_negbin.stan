// Cairns-Blake-Dowd sous loi binomiale negative.
// Les deux facteurs temporels suivent une marche aleatoire bivariee.
// La pente est echantillonnee par ecart-type d'age pour equilibrer
// numeriquement le niveau et la pente. La matrice kappa restitue ensuite
// le coefficient par annee d'age de la formule CBD usuelle.

data {
  int<lower=2> A;
  int<lower=2> T;
  array[A, T] int<lower=0> D;
  matrix[A, T] log_E;
  vector[A] age_centered;
  real<lower=0> age_scale;
}

transformed data {
  int N = A * T;
  array[N] int D_long;
  vector[N] log_E_long;
  for (t in 1:T) {
    for (a in 1:A) {
      int i = (t - 1) * A + a;
      D_long[i] = D[a, t];
      log_E_long[i] = log_E[a, t];
    }
  }
}

parameters {
  matrix[2, T] kappa_scaled;
  vector[2] drift_scaled;
  vector<lower=0>[2] sigma_kappa_scaled;
  cholesky_factor_corr[2] L_Omega;
  real<lower=0> inv_phi;
}

transformed parameters {
  matrix[2, T] kappa;
  vector[2] drift;
  vector<lower=0>[2] sigma_kappa;
  matrix[2, 2] L_Sigma =
    diag_pre_multiply(sigma_kappa_scaled, L_Omega);
  corr_matrix[2] Omega = multiply_lower_tri_self_transpose(L_Omega);
  real<lower=0> phi = inv(inv_phi);

  drift[1] = drift_scaled[1];
  drift[2] = drift_scaled[2] / age_scale;
  sigma_kappa[1] = sigma_kappa_scaled[1];
  sigma_kappa[2] = sigma_kappa_scaled[2] / age_scale;

  // CBD n'a pas d'intercept alpha_x : ses etats initiaux sont estimes.
  for (t in 1:T) {
    kappa[1, t] = kappa_scaled[1, t];
    kappa[2, t] = kappa_scaled[2, t] / age_scale;
  }
}

model {
  vector[A * T] eta;

  // Priors faibles mais places sur l'echelle des log-taux de mortalite.
  kappa_scaled[1, 1] ~ normal(-4, 2);
  kappa_scaled[2, 1] ~ normal(0, 2);
  drift_scaled ~ normal(0, 0.1);
  sigma_kappa_scaled ~ exponential(10);
  L_Omega ~ lkj_corr_cholesky(2);
  // Forme centree : chaque etat est fortement informe par 41 ages.
  for (t in 2:T) {
    kappa_scaled[, t] ~ multi_normal_cholesky(
      kappa_scaled[, t - 1] + drift_scaled,
      L_Sigma
    );
  }
  inv_phi ~ normal(0, 1);

  for (t in 1:T) {
    for (a in 1:A) {
      int i = (t - 1) * A + a;
      eta[i] = kappa[1, t] + age_centered[a] * kappa[2, t];
    }
  }
  D_long ~ neg_binomial_2_log(log_E_long + eta, phi);
}
