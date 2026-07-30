// Lee-Carter sous loi binomiale negative.
// Identification : kappa[1] = 0 et beta appartient au simplexe.

data {
  int<lower=2> A;
  int<lower=2> T;
  array[A, T] int<lower=0> D;
  matrix[A, T] log_E;
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
  vector[A] alpha;
  simplex[A] beta;
  real drift;
  real<lower=0> sigma_kappa;
  vector[T - 1] kappa_free;
  real<lower=0> inv_phi;
}

transformed parameters {
  vector[T] kappa;
  real<lower=0> phi = inv(inv_phi);

  kappa[1] = 0;
  kappa[2:T] = kappa_free;
}

model {
  vector[A * T] eta;

  alpha ~ normal(-5, 3);
  beta ~ dirichlet(rep_vector(1.0, A));
  drift ~ normal(0, 2);
  sigma_kappa ~ exponential(0.5);
  for (t in 2:T) {
    kappa[t] ~ normal(kappa[t - 1] + drift, sigma_kappa);
  }
  inv_phi ~ normal(0, 1);

  for (t in 1:T) {
    for (a in 1:A) {
      int i = (t - 1) * A + a;
      eta[i] = alpha[a] + beta[a] * kappa[t];
    }
  }
  D_long ~ neg_binomial_2_log(log_E_long + eta, phi);
}
