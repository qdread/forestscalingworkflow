# Script to fit to binned richness values, first for all trees then with a mixed model with coefficients for each FG
# Using both diameter and light/area

# The following procedure is done in this script:
# 1. Load data 
# 2. Load and fit models
# 3. Extract output
#   - parameter values
#   - fitted values
#   - Bayesian R-squared
# 4. Create clean summary tables
# 5. Create plotting data CSVs


# Load data ----------------------------------------------------

library(tidyverse)
library(forestscaling)
library(rstan)
options(mc.cores = 3)

load('data/rawdataobj_alternativecluster.r')

# Load 1995 bin data (to make sure we're using the same bin breaks as other figs)
load('data/data_binned/bin_object_singleyear.RData')

# Create binned richness data

# Use existing bin bounds
bin_bounds <- fastslow_stats_bydiam_byyear[1:20, c('bin_midpoint', 'bin_min', 'bin_max')]
bin_bounds_light <- fastslow_stats_bylight_byyear[1:20, c('bin_midpoint', 'bin_min', 'bin_max')]
# Get richness of each FG by each bin

bin_x_fg <- expand_grid(fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_x_fg_light <- expand_grid(fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds_light %>% mutate(bin = 1:20))

# 1995 data (132,982 trees)
dat <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

# 1995 data with light values (113,651 trees)
dat_light <- dat %>%
  filter(!is.na(light_received_byarea))

bin_x_fg <- bin_x_fg %>%
  cbind(pmap_dfr(bin_x_fg, function(fg, bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$fg %in% fg & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

bin_x_fg_light <- bin_x_fg_light %>%
  cbind(pmap_dfr(bin_x_fg_light, function(fg, bin_min, bin_max, ...) {
    sp_ids <- as.character(dat_light$sp[dat_light$fg %in% fg & dat_light$light_received_byarea >= bin_min & dat_light$light_received_byarea < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Set up binned values for "all" 

bin_all <- data.frame(fg = 'all', bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_all_light <- data.frame(fg = 'all', bin = 1:20) %>%
  left_join(bin_bounds_light %>% mutate(bin = 1:20))

bin_all <- bin_all %>%
  cbind(pmap_dfr(bin_all, function(bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

bin_all_light <- bin_all_light %>%
  cbind(pmap_dfr(bin_all_light, function(bin_min, bin_max, ...) {
    sp_ids <- as.character(dat_light$sp[dat_light$light_received_byarea >= bin_min & dat_light$light_received_byarea < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Create model input data for all trees
dat_all <- with(bin_all %>% filter(richness > 0), list(x = bin_midpoint, y = richness_by_bin_width, N = nrow(bin_all)))
dat_all_light <- with(bin_all_light %>% filter(richness > 0), list(x = bin_midpoint, y = richness_by_bin_width, N = nrow(bin_all_light)))

# Create model input data for each FG
bin_x_fg_use <- bin_x_fg %>% filter(richness > 0, !fg %in% 'unclassified')
dat_fg <- with(bin_x_fg_use, list(x = bin_midpoint, y = richness_by_bin_width, fg = as.numeric(factor(fg)), N = nrow(bin_x_fg_use), M = 5))
bin_x_fg_light_use <- bin_x_fg_light %>% filter(richness > 0, !fg %in% 'unclassified')
dat_fg_light <- with(bin_x_fg_light_use, list(x = bin_midpoint, y = richness_by_bin_width, fg = as.numeric(factor(fg)), N = nrow(bin_x_fg_light_use), M = 5))

# Load and fit models -----------------------------------------------------

mod_2seg_linear <- stan_model('model_scripts/richness_2segment_linearmodel.stan')
mod_3seg_linear <- stan_model('model_scripts/richness_3segment_linearmodel.stan')

mod_2seg_mixed <- stan_model('model_scripts/richness_2segment_mixedmodel.stan')
mod_3seg_mixed <- stan_model('model_scripts/richness_3segment_mixedmodel.stan')

# Note: warnings about Rhat = NA are an artifact of how the indicator variable is defined. You can ignore them.

fit_2seg_linear_all <- sampling(mod_2seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 222)
fit_3seg_linear_all <- sampling(mod_3seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 111)

fit_2seg_linear_all_light <- sampling(mod_2seg_linear, data = dat_all_light, chains = 3, iter = 5000, warmup = 4000, seed = 222)
fit_3seg_linear_all_light <- sampling(mod_3seg_linear, data = dat_all_light, chains = 3, iter = 5000, warmup = 4000, seed = 111)

fit_2seg_mixed_all <- sampling(mod_2seg_mixed, data = dat_fg, chains = 3, iter = 5000, warmup = 4000, seed = 333)
# More iterations needed to converge 3 segment mixed model, and increased treedepth. Works OK.
fit_3seg_mixed_all <- sampling(mod_3seg_mixed, data = dat_fg, chains = 3, iter = 10000, warmup = 9000, seed = 4444, control = list(max_treedepth = 20))

fit_2seg_mixed_all_light <- sampling(mod_2seg_mixed, data = dat_fg_light, chains = 3, iter = 5000, warmup = 4000, seed = 555)
fit_3seg_mixed_all_light <- sampling(mod_3seg_mixed, data = dat_fg_light, chains = 3, iter = 10000, warmup = 9000, seed = 666, control = list(max_treedepth = 20))

# Extract output ----------------------------------------------------------

# Functions for fitted values
# ===========================

twoseg_log <- function(x, alpha, beta, tau) {
  x2 <- ifelse(log10(x) < tau, 0, 1)
  10 ^ (alpha + beta[1] * log10(x) + beta[2] * (log10(x) - tau) * x2)
}

threeseg_log <- function(x, alpha, beta, tau_low, tau_high) {
  x2a <- ifelse(log10(x) < tau_low, 0, 1)
  x2b <- ifelse(log10(x) < tau_high, 0, 1)
  10 ^ (alpha + beta[1] * log10(x) + beta[2] * (log10(x) - tau_low) * x2a + beta[3] * (log10(x) - tau_high) * x2b)
}

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))
light_pred <- exp(seq(log(1), log(412), length.out = 101))
qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 

# Only extract from three segment model.
# Extract from the all-tree linear model, and the by-group mixed model, for both diameter and light/area.

# Parameter values
# ================

# Create vector with parameter names to extract
# In linear model, beta is a 3 element vector.
# Only get the group level coefficients from the mixed model.
par_names_3seg_linear <- c('alpha', 'beta', 'tau_low', 'tau_high')
par_names_3seg_mixed <- c('coef_alpha', 'coef_beta_low', 'coef_beta_mid', 'coef_beta_high', 'coef_tau_low', 'coef_tau_high')
df_col_names <- c('fg', 'parameter', 'mean', 'se_mean', 'sd', 'q025', 'q25', 'q50', 'q75', 'q975', 'n_eff', 'Rhat')

params_3seg_diam_all <- data.frame(fg = 'all',
                                   parameter = c('alpha','beta_low','beta_mid','beta_high','tau_low','tau_high'),
                                   summary(fit_3seg_linear_all, pars = par_names_3seg_linear)[[1]]) %>%
  setNames(df_col_names)
params_3seg_diam_mixed <- data.frame(fg = 1:5,
                                     parameter = rep(c('alpha','beta_low','beta_mid','beta_high','tau_low','tau_high'), each = 5),
                                     summary(fit_3seg_mixed_all, pars = par_names_3seg_mixed)[[1]]) %>%
  setNames(df_col_names)
params_3seg_diam <- rbind(params_3seg_diam_all, params_3seg_diam_mixed)

params_3seg_light_all <- data.frame(fg = 'all',
                                    parameter = c('alpha','beta_low','beta_mid','beta_high','tau_low','tau_high'),
                                   summary(fit_3seg_linear_all_light, pars = par_names_3seg_linear)[[1]]) %>%
  setNames(df_col_names)
params_3seg_light_mixed <- data.frame(fg = 1:5,
                                     parameter = rep(c('alpha','beta_low','beta_mid','beta_high','tau_low','tau_high'), each = 5),
                                     summary(fit_3seg_mixed_all_light, pars = par_names_3seg_mixed)[[1]]) %>%
  setNames(df_col_names)
params_3seg_light <- rbind(params_3seg_light_all, params_3seg_light_mixed)

# Fitted and predicted values
# ===========================

# Extract needed parameter estimates from each iteration, to be used for fitted values, fitted slopes, R2, and corr. factor.
par_draws_linear_diam <- extract(fit_3seg_linear_all, pars = par_names_3seg_linear)
par_draws_linear_light <- extract(fit_3seg_linear_all_light, pars = par_names_3seg_linear)
par_draws_mixed_diam <- extract(fit_3seg_mixed_all, pars = par_names_3seg_mixed)
par_draws_mixed_light <- extract(fit_3seg_mixed_all_light, pars = par_names_3seg_mixed)

# Reshape the parameter draws into data frames
par_draws_linear_diam <- do.call(cbind, par_draws_linear_diam) %>%
  as.data.frame %>% setNames(c('alpha', 'beta_low', 'beta_mid', 'beta_high', 'tau_low', 'tau_high'))
par_draws_linear_light <- do.call(cbind, par_draws_linear_light) %>%
  as.data.frame %>% setNames(c('alpha', 'beta_low', 'beta_mid', 'beta_high', 'tau_low', 'tau_high'))

par_draws_mixed_diam <- imap_dfr(par_draws_mixed_diam, ~ data.frame(iter = 1:nrow(.x), parameter = .y, .x)) %>%
  pivot_longer(-c(iter, parameter), names_to = 'fg') %>%
  pivot_wider(id_cols = c(iter, fg), names_from = parameter, values_from = value) %>%
  mutate(fg = as.numeric(substr(fg, 2, 2)))
par_draws_mixed_light <- imap_dfr(par_draws_mixed_light, ~ data.frame(iter = 1:nrow(.x), parameter = .y, .x)) %>%
  pivot_longer(-c(iter, parameter), names_to = 'fg') %>%
  pivot_wider(id_cols = c(iter, fg), names_from = parameter, values_from = value) %>%
  mutate(fg = as.numeric(substr(fg, 2, 2)))

# Get fitted values for each iteration (and FG)
fitted_linear_diam <- par_draws_linear_diam %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = dbh_pred, y_fit = threeseg_log(x = dbh_pred, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))
fitted_linear_light <- par_draws_linear_light %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = light_pred, y_fit = threeseg_log(x = light_pred, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))


fitted_mixed_diam <- par_draws_mixed_diam %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = dbh_pred, y_fit = threeseg_log(x = dbh_pred, alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)))
fitted_mixed_light <- par_draws_mixed_light %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = light_pred, y_fit = threeseg_log(x = light_pred, alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)))


# Get quantiles of the fitted values
df_col_names <- c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')

fittedquant_linear_diam <- fitted_linear_diam %>%
  group_by(x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('dbh', df_col_names))

fittedquant_linear_light <- fitted_linear_light %>%
  group_by(x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('lightarea', df_col_names))

fittedquant_mixed_diam <- fitted_mixed_diam %>%
  group_by(fg, x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('fg', 'dbh', df_col_names))

fittedquant_mixed_light <- fitted_mixed_light %>%
  group_by(fg, x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('fg', 'lightarea', df_col_names))

# correction factor for log-log regression plot
# =============================================

# Get fitted values for every data point at every iteration (posterior estimates of linear predictor)
linpred_linear_diam <- par_draws_linear_diam %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = dat_all$x, y = dat_all$y, y_fit = threeseg_log(x = dat_all$x, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))
linpred_linear_light <- par_draws_linear_light %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = dat_all_light$x, y = dat_all_light$y, y_fit = threeseg_log(x = dat_all_light$x, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))

linpred_mixed_diam <- par_draws_mixed_diam %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = dat_fg$x[dat_fg$fg == .$fg[1]],
                            y = dat_fg$y[dat_fg$fg == .$fg[1]],
                            y_fit = threeseg_log(x =  dat_fg$x[dat_fg$fg == .$fg[1]], alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)),
               keep = TRUE)
linpred_mixed_light <- par_draws_mixed_light %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = dat_fg_light$x[dat_fg_light$fg == .$fg[1]],
                            y = dat_fg_light$y[dat_fg_light$fg == .$fg[1]],
                            y_fit = threeseg_log(x =  dat_fg_light$x[dat_fg_light$fg == .$fg[1]], alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)),
               keep = TRUE)

# Take log of fitted values and get residuals by subtracting log of observed y from log of fitted y
linpred_linear_diam <- linpred_linear_diam %>%
  mutate(resid = log(y_fit) - log(y))
linpred_linear_light <- linpred_linear_light %>%
  mutate(resid = log(y_fit) - log(y))

linpred_mixed_diam <- linpred_mixed_diam %>%
  mutate(resid = log(y_fit) - log(y))
linpred_mixed_light <- linpred_mixed_light %>%
  mutate(resid = log(y_fit) - log(y))

# Bias correction factor: take sum of squared residuals
# Get standard error of estimates, using number of parameters to correct the degrees of freedom
# Use formula to get correction factor

k_linearmodel <- 7
k_mixedmodel <- 13

cf_linear_diam <- linpred_linear_diam %>%
  group_by(iter) %>%
  summarize(ssq_resid = sum(resid^2),
            n = n()) %>%
  mutate(see = (ssq_resid / (n - k_linearmodel))^0.5,
         cf = exp((see^2)/2))
cf_linear_light <- linpred_linear_light %>%
  group_by(iter) %>%
  summarize(ssq_resid = sum(resid^2),
            n = n()) %>%
  mutate(see = (ssq_resid / (n - k_linearmodel))^0.5,
         cf = exp((see^2)/2))

cf_mixed_diam <- linpred_mixed_diam %>%
  group_by(iter) %>%
  summarize(ssq_resid = sum(resid^2),
            n = n()) %>%
  mutate(see = (ssq_resid / (n - k_mixedmodel))^0.5,
         cf = exp((see^2)/2))
cf_mixed_light <- linpred_mixed_light %>%
  group_by(iter) %>%
  summarize(ssq_resid = sum(resid^2),
            n = n()) %>%
  mutate(see = (ssq_resid / (n - k_mixedmodel))^0.5,
         cf = exp((see^2)/2))
  
# Extract quantiles of correction factor


# Write the extracted output to files
# ===================================

# Create clean summary tables ---------------------------------------------

# FIXME put the code here

# Create plotting data ----------------------------------------------------

# FIXME put the code here
