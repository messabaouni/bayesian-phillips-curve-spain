# Regional Phillips Curve in Spain: A Hierarchical Bayesian Approach

Bayesian Data Analysis project (MESIO UPC-UB, 2025), written in R and Stan.
It tests whether the Phillips curve holds across Spain's autonomous communities, by relating inflation to unemployment (`paro`), GDP, GDP growth and unit labour cost (`salarycost`) with regional effects.

## Models (Stan)
1. `modelo_simple.stan`: Normal likelihood, region-varying intercept and common slopes
2. `modelo_intermedio.stan`: also lets the unemployment slope vary by region
3. `modelo_gamma.stan`: Gamma likelihood with a log link for positive, right-skewed inflation, hierarchical slopes, a quadratic unemployment term and an unemployment × GDP interaction

The analysis includes convergence traces, posterior distributions by region, and posterior predictive checks (observed vs predicted, residuals, MAE and RMSE per region).

**Main findings:** economic growth is clearly associated with higher inflation, while the unemployment effect is weaker and varies by region. The Phillips curve does not hold uniformly across regions.

## Reading the report
The report builds the models in steps, from a simple model to the full hierarchical Gamma model, so that each added feature (region-specific slopes, non-linear terms, a Gamma likelihood) can be justified. Diagnostics use trace plots and posterior predictive checks. The full write-up is in `report/report.pdf` (in Spanish).

## Data note
The regional data (GDP, labour costs, inflation and unemployment by autonomous community, 2018 onwards) comes from Spanish public statistics. The exact source tables are not documented in the repository, so check the original providers' terms before reusing the data.

## Structure
```
data/data.csv    regional panel (region, year, GDP, labour cost, inflation, unemployment, GDP growth)
stan/            the three Stan models
R/analysis.R     model fitting script
report/report.Rmd  full report source (in Spanish)
report/report.pdf  rendered report
```

## How to run
From the repository root:
```r
source("R/analysis.R")
```
Render the report with `rmarkdown::render("report/report.Rmd")`.

Packages: rstan, bayesplot, posterior, dplyr, tidyr, ggplot2, patchwork, readr, knitr, tibble.

## Authors
Julia Anglada and Mohammed Essabaouni
