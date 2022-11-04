# Script for additional mortality versus diameter fit
# J-curve functional form
# Each nonlinear parameter has functional group as fixed effect
# 17 Sep 2022

# Modified 04 November 2022: Use consistent method for calculating mortality.

# Setup -------------------------------------------------------------------

library(dplyr)
library(readr)
library(brms)
library(tidybayes)
library(ggplot2)
library(tidyr)
library(purrr)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'never')

load('data/rawdataobj_alternativecluster.r')

# Load 1995 bin data (to make sure we're using the same bin breaks as other figs)
load('data/data_binned/bin_object_singleyear.RData')

# Generate 1990->1995 mortality.
dat_mort <- map2_dfr(alltreedat, seq(1985, 2010, by = 5), function(x,y) {
  x %>%
    select(treeID, sp, fg, dbh_corr) %>%
    mutate(year = y)
}) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  group_by(treeID) %>%
  mutate(last_year = max(year)) %>% 
  ungroup %>%
  mutate(died = year == last_year)

mort_data <- dat_mort %>%
  filter(year == 1990, !fg %in% 'unclassified') %>%
  select(fg, died, dbh = dbh_corr)


# Fit model ---------------------------------------------------------------

mort_jcurve_fixef_fit <- brm(
  bf(
    died ~ alpha + -exp(beta) * log(dbh) * exp(gamma * log(dbh)),
    alpha ~ 0 + fg,
    beta ~ 0 + fg,
    gamma ~ 0 + fg,
    nl = TRUE
  ),
  data = mort_data, family = bernoulli(link = 'logit'),
  prior = c(
    prior(normal(0, 2), nlpar = alpha),
    prior(normal(0, 2), nlpar = beta),
    prior(normal(0, 2), nlpar = gamma)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 27701,
  file = '~/temp/forestlight/mort_jcurve_fixef_brmfit_new'
)


# Generate predicted values -----------------------------------------------

# Read pre-created mortality bins for plotting to compare to fitted values
binned_data <- read_csv('data/data_binned/additional_bins_fg_year.csv')

# Prediction grid: dbh x fg
# Limit ranges to the observed data points
mort_max <- binned_data %>%
  filter(abundance > 20, fg %in% paste0('fg', 1:5)) %>%
  group_by(fg, year) %>%
  filter(bin_midpoint == max(bin_midpoint))

pred_dat <- mort_max %>%
  filter(year == 1995) %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dbh = exp(seq(log(1), log(.$bin_midpoint), length.out = 101))))

# Get expected values of posterior means (epred)
jcurve_pred <- pred_dat %>%
  add_epred_draws(mort_jcurve_fixef_fit) %>%
  mutate(.epred = -log(1-.epred)/5)

# Generate parameter tables -----------------------------------------------

mort_jcurve_params <- mort_jcurve_fixef_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  separate(.variable, into = c('b', 'parameter', 'fg')) %>%
  select(-b) %>%
  mutate(fg = as.integer(substr(fg, 5, 5))) %>%
  mutate(.value = if_else(parameter == 'beta', -exp(.value), exp(.value)))

mort_jcurve_quant <- mort_jcurve_params %>%
  group_by(parameter, fg) %>%
  median_qi(.width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(parameter, fg, .value), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, parameter, .value, .lower_0.95, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'parameter', 'median', 'q025', 'q05', 'q25', 'q75', 'q95', 'q975')) %>%
  mutate(fg = rep(c('fast', 'large pioneer', 'slow', 'small breeder', 'medium'), 3),
         parameter = rep(c('alpha (intercept)', 'beta (slope at small size)', 'gamma (parameter determining J-curve location)'), each = 5))

# Generate plotting data --------------------------------------------------

mort_jcurve_plotting_data <- jcurve_pred %>%
  group_by(fg, dbh) %>%
  median_qi(.epred, .width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(fg, dbh, .epred), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, dbh, .lower_0.95, .lower_0.9, .lower_0.5, .epred, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')) 

write_csv(mort_jcurve_quant, 'data/clean_summary_tables/clean_mortalitybydiameter_parameters.csv')
write_csv(mort_jcurve_plotting_data, 'data/data_forplotting/fitted_mortalitybydiameter.csv')
