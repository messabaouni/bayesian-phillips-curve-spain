data {
  int<lower=1> N;
  int<lower=1> R;
  int<lower=1, upper=R> region[N];
  vector[N] paro;
  vector[N] pib;
  vector[N] pib_growth;
  vector[N] salarycost;
  vector[N] inflation;
}

parameters {
  real beta_paro;
  real beta_pib;
  real beta_growth;
  real beta_salary;

  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[R] alpha_raw;

  real<lower=0> sigma;
}

transformed parameters {
  vector[R] alpha;
  alpha = mu_alpha + sigma_alpha * alpha_raw;
}

model {
  // Priors
  beta_paro ~ normal(0, 1);
  beta_pib ~ normal(0, 1);
  beta_growth ~ normal(0, 1);
  beta_salary ~ normal(0, 1);

  mu_alpha ~ normal(log(1.02), 1);
  sigma_alpha ~ normal(0, 1);
  alpha_raw ~ normal(0, 1);

  sigma ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    real mu;
    mu = alpha[region[n]]
        + beta_paro * paro[n]
        + beta_pib * pib[n]
        + beta_growth * pib_growth[n]
        + beta_salary * salarycost[n];

    inflation[n] ~ normal(mu, sigma);
  }
}

