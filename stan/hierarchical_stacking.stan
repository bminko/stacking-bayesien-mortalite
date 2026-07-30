// Meta-modele de stacking contextuel hierarchique.
// X contient, dans cet ordre : age, age^2 centre, horizon, interaction.

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;
  matrix[N, K] log_p;
  matrix[N, P] X;
  vector<lower=0>[N] obs_weight;
  vector<lower=0>[P] tau_scale;
}

parameters {
  vector[K - 1] alpha;
  matrix[K - 1, P] z_beta;
  vector<lower=0>[P] tau;
}

transformed parameters {
  matrix[K - 1, P] beta;
  for (p in 1:P) {
    beta[, p] = tau[p] * z_beta[, p];
  }
}

model {
  alpha ~ normal(0, 1.5);
  to_vector(z_beta) ~ std_normal();
  tau ~ normal(0, tau_scale);

  for (n in 1:N) {
    vector[K] score;
    vector[K] log_weight;
    score[1:(K - 1)] = alpha + beta * to_vector(X[n]');
    score[K] = 0;
    log_weight = score - log_sum_exp(score);
    target += obs_weight[n] * log_sum_exp(log_weight + to_vector(log_p[n]'));
  }
}
