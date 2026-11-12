
data {
  int<lower=1> N;
  int<lower=1> R;
  int<lower=1, upper=R> region[N];
  vector[N] paro;
  vector[N] pib;
  vector[N] pib_growth;
  vector[N] salarycost;
  vector<lower=0>[N] inflation;
}

parameters {
  // Hiperparámetros globales
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real mu_beta_paro;
  real mu_beta_pib;
  real mu_beta_growth;
  real mu_beta_salary;

  real<lower=0> sigma_beta_paro;
  real<lower=0> sigma_beta_pib;
  real<lower=0> sigma_beta_growth;
  real<lower=0> sigma_beta_salary;

  // Efectos por región
  vector[R] alpha;
  vector[R] beta_paro;
  vector[R] beta_pib;
  vector[R] beta_growth;
  vector[R] beta_salary;

  // Parámetros funcionales adicionales
  real beta_paro_sq;
  real beta_interaction_paro_pib;

  // Parámetro de forma de la distribución Gamma
  real<lower=0> shape;
}

transformed parameters {
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = alpha[region[n]]
          + beta_paro[region[n]] * paro[n]
          + beta_pib[region[n]] * pib[n]
          + beta_growth[region[n]] * pib_growth[n]
          + beta_salary[region[n]] * salarycost[n]
          + beta_paro_sq * square(paro[n])
          + beta_interaction_paro_pib * paro[n] * pib[n];
    mu[n] = exp(mu[n]);
  }
}

model {
  // Priors hiperparámetros
  mu_alpha ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);

  mu_beta_paro ~ normal(0, 1);
  mu_beta_pib ~ normal(0, 1);
  mu_beta_growth ~ normal(0, 1);
  mu_beta_salary ~ normal(0, 1);

  sigma_beta_paro ~ normal(0, 1);
  sigma_beta_pib ~ normal(0, 1);
  sigma_beta_growth ~ normal(0, 1);
  sigma_beta_salary ~ normal(0, 1);

  // Priors por región
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta_paro ~ normal(mu_beta_paro, sigma_beta_paro);
  beta_pib ~ normal(mu_beta_pib, sigma_beta_pib);
  beta_growth ~ normal(mu_beta_growth, sigma_beta_growth);
  beta_salary ~ normal(mu_beta_salary, sigma_beta_salary);

  // Priors funcionales
  beta_paro_sq ~ normal(0, 1);
  beta_interaction_paro_pib ~ normal(0, 1);

  shape ~ gamma(2, 1); // prior para forma

  // Likelihood Gamma (usando parametrización shape y mean)
  for (n in 1:N) {
    inflation[n] ~ gamma(shape, shape / mu[n]);
  }
}
