
suppressPackageStartupMessages({
  library(cmdstanr)
})


mod <- cmdstan_model(stan_file = "NonStrat.stan")
mod <- cmdstan_model(stan_file = "Strat_noKappa.stan")
mod <- cmdstan_model(stan_file = "Strat_oneKappa.stan")