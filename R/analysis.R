
# paquetes
library(readr)
library(dplyr)
library(ggplot2)
library(knitr)
library(tidyr)    
library(bayesplot) 
library(rstan)
library(bayesplot)
library(patchwork)
library(posterior)
library(tibble)

# Lectura de datos
data <- read.csv("data/data.csv", sep = ";")
data_std <- data %>%
  mutate(region_id = as.integer(as.factor(region))) %>%
  mutate(
    paro = scale(paro)[,1],
    pib = scale(pib)[,1],
    pib_growth = scale(pib_growth)[,1],
    salarycost = scale(salarycost)[,1]
  ) %>%
  filter(inflation > 0)

stan_data <- list(
  N = nrow(data_std),
  R = length(unique(data_std$region_id)),
  region = data_std$region_id,
  paro = data_std$paro,
  pib = data_std$pib,
  pib_growth = data_std$pib_growth,
  salarycost = data_std$salarycost,
  inflation = data_std$inflation
)

# Modelos en Stan
simple_fit <- stan(
  file = "stan/modelo_simple.stan",
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

intermedio_fit <- stan(
  file = "stan/modelo_intermedio.stan",
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

gamma_fit <- stan(
  file = "stan/modelo_gamma.stan",
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# Trazas
theme_trace <- theme_minimal(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"))

trace_simple <- mcmc_trace(simple_fit,
                           pars = c("beta_paro", "beta_pib", "beta_growth", "beta_salary")) +
  ggtitle("Modelo simple") +
  theme_trace

trace_intermedio <- mcmc_trace(intermedio_fit,
                                pars = c("mu_beta_paro", "beta_pib", "beta_growth", "beta_salary")) +
  ggtitle("Modelo intermedio") +
  theme_trace

trace_gamma <- mcmc_trace(gamma_fit,
                          pars = c("mu_beta_paro", "mu_beta_pib", "mu_beta_growth", "mu_beta_salary")) +
  ggtitle("Modelo gamma") +
  theme_trace

(trace_simple / trace_intermedio / trace_gamma) +
  plot_annotation(title = "Trazas de convergencia comparadas entre modelos")

# Estimaciones
summary_table <- summary(gamma_fit, pars = c("mu_alpha", "mu_beta_paro", "mu_beta_pib", "mu_beta_growth", "mu_beta_salary",
                          "shape", "beta_paro_sq", "beta_interaction_paro_pib"))$summary
print(round(summary_table, 3))

plot(gamma_fit, pars = c("mu_alpha", "mu_beta_paro", "mu_beta_pib", "mu_beta_growth", "mu_beta_salary",
                         "beta_paro_sq", "beta_interaction_paro_pib"))

# Intervalos globales
posterior_df <- as_draws_df(gamma_fit)

mcmc_areas(
  posterior_df,
  pars = c("mu_alpha", "mu_beta_paro", "mu_beta_pib", "mu_beta_growth", "mu_beta_salary"),
  prob = 0.95
) +
  ggtitle("Intervalos de credibilidad del 95% para parámetros globales") +
  theme_minimal()

# Interceptos regionales
alpha_draws <- posterior_df %>%
  select(starts_with("alpha[")) %>%
  pivot_longer(cols = everything(), names_to = "region", values_to = "alpha")

ggplot(alpha_draws, aes(x = region, y = alpha)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Distribución posterior de los interceptos regionales",
       x = "Región", y = expression(alpha[r])) +
  theme_minimal()

# Efectos regionales del crecimiento del PIB
g1 <- mcmc_areas(
  posterior_df,
  pars = paste0("beta_growth[", 1:10, "]"),
  prob = 0.95
) +
  ggtitle("Regiones 1–10") +
  theme_minimal()

g2 <- mcmc_areas(
  posterior_df,
  pars = paste0("beta_growth[", 11:19, "]"),
  prob = 0.95
) +
  ggtitle("Regiones 11–19") +
  theme_minimal()

(g1 | g2) +
  plot_annotation(
    title = "Distribución posterior de los efectos regionales del crecimiento del PIB"
  )

# Efectos regionales del paro
mcmc_areas(
  posterior_df,
  pars = paste0("beta_paro[", 1:19, "]"),
  prob = 0.95
) +
  ggtitle("Distribución posterior de los efectos regionales del paro") +
  theme_minimal()

# Posterior predictiva
posterior_df <- as.data.frame(gamma_fit)
n_sims <- 5000
N <- stan_data$N
region <- stan_data$region

inflation_rep <- matrix(NA, nrow = n_sims, ncol = N)

for (s in 1:n_sims) {
  alpha     <- posterior_df[s, grep("^alpha\[", names(posterior_df))]
  beta_paro <- posterior_df$mu_beta_paro[s]
  beta_pib  <- posterior_df$mu_beta_pib[s]
  beta_grow <- posterior_df$mu_beta_growth[s]
  beta_sal  <- posterior_df$mu_beta_salary[s]
  shape     <- posterior_df$shape[s]

  mu <- as.numeric(alpha[region]) +
        beta_paro * stan_data$paro +
        beta_pib  * stan_data$pib +
        beta_grow * stan_data$pib_growth +
        beta_sal  * stan_data$salarycost

  mu <- exp(mu)
  rate <- shape / mu
  inflation_rep[s, ] <- rgamma(N, shape = shape, rate = rate)
}

pred_mean <- apply(inflation_rep, 2, mean)
lower_95  <- apply(inflation_rep, 2, quantile, 0.025)
upper_95  <- apply(inflation_rep, 2, quantile, 0.975)

ppc_data <- tibble::tibble(
  y_obs = stan_data$inflation,
  y_pred = pred_mean,
  lower = lower_95,
  upper = upper_95
)

# Gráfico observaciones vs predicciones
ppc_data <- ppc_data %>% filter(y_pred < 15)

ggplot(ppc_data, aes(x = y_obs, y = y_pred)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior predictiva: Observado vs Predicho (sin valores extremos)",
    x = "Inflación observada",
    y = "Media predicha con IC95%"
  ) +
  theme_minimal()

# Gráfico residuos
ppc_data <- ppc_data %>%
  mutate(residual = y_obs - y_pred)

ggplot(ppc_data, aes(x = y_pred, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    title = "Residuos vs Predicción",
    x = "Inflación predicha",
    y = "Residual (observado - predicho)"
  ) +
  theme_minimal()

# MAE y RMSE por región
ppc_df <- tibble(
  region = stan_data$region,
  y_obs = stan_data$inflation,
  y_pred = pred_mean
)

errores_region <- ppc_df %>%
  group_by(region) %>%
  summarise(
    MAE = mean(abs(y_obs - y_pred)),
    RMSE = sqrt(mean((y_obs - y_pred)^2))
  ) %>%
  arrange(region)

print(errores_region)
