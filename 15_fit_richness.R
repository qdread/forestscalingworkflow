# Script to fit to binned richness values, first for all trees then with a mixed model with coefficients for each FG
# Using both diameter and light/area

# The following procedure is done in this script:
# 1. Load data 
# 2. Load and fit models
# 3. Extract output
#   - parameter values
#   - fitted values
#   - correction factor for log log regression plots
# 4. Create clean summary tables
# 5. Create plotting data CSVs
# 6. Create fitted values and fitted slopes for the ratio plots as well, output them to CSVs (one for plotting and one clean)
# 7. Binned abundance versus binned richness values, fit a regression line to them (mixed-effects model)


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
# Get richness of each FG by each bin

bin_x_fg <- expand_grid(fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

# 1995 data (132,982 trees)
dat <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

bin_x_fg <- bin_x_fg %>%
  cbind(pmap_dfr(bin_x_fg, function(fg, bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$fg %in% fg & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Set up binned values for "all" 

bin_all <- data.frame(fg = 'all', bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_all <- bin_all %>%
  cbind(pmap_dfr(bin_all, function(bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Create model input data for all trees
dat_all <- with(bin_all %>% filter(richness > 0), list(x = bin_midpoint, y = richness_by_bin_width, N = nrow(bin_all)))

# Create model input data for each FG
bin_x_fg_use <- bin_x_fg %>% filter(richness > 0, !fg %in% 'unclassified')
dat_fg <- with(bin_x_fg_use, list(x = bin_midpoint, y = richness_by_bin_width, fg = as.numeric(factor(fg)), N = nrow(bin_x_fg_use), M = 5))

# Load and fit models -----------------------------------------------------

mod_2seg_linear <- stan_model('model_scripts/richness_2segment_linearmodel.stan')
mod_3seg_linear <- stan_model('model_scripts/richness_3segment_linearmodel.stan')

mod_2seg_mixed <- stan_model('model_scripts/richness_2segment_mixedmodel.stan')
mod_3seg_mixed <- stan_model('model_scripts/richness_3segment_mixedmodel.stan')

# Note: warnings about Rhat = NA are an artifact of how the indicator variable is defined. You can ignore them.

fit_2seg_linear_all <- sampling(mod_2seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 222)
fit_3seg_linear_all <- sampling(mod_3seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 111)

fit_2seg_mixed_all <- sampling(mod_2seg_mixed, data = dat_fg, chains = 3, iter = 5000, warmup = 4000, seed = 333)
# More iterations needed to converge 3 segment mixed model, and increased treedepth. Works OK.
fit_3seg_mixed_all <- sampling(mod_3seg_mixed, data = dat_fg, chains = 3, iter = 10000, warmup = 9000, seed = 4444, control = list(max_treedepth = 20))

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
qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 

# Only extract from three segment model.
# Extract from the all-tree linear model, and the by-group mixed model

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

# Fitted and predicted values
# ===========================

# Extract needed parameter estimates from each iteration, to be used for fitted values, fitted slopes, R2, and corr. factor.
par_draws_linear_diam <- extract(fit_3seg_linear_all, pars = par_names_3seg_linear)
par_draws_mixed_diam <- extract(fit_3seg_mixed_all, pars = par_names_3seg_mixed)

# Reshape the parameter draws into data frames
par_draws_linear_diam <- do.call(cbind, par_draws_linear_diam) %>%
  as.data.frame %>% setNames(c('alpha', 'beta_low', 'beta_mid', 'beta_high', 'tau_low', 'tau_high'))

par_draws_mixed_diam <- imap_dfr(par_draws_mixed_diam, ~ data.frame(iter = 1:nrow(.x), parameter = .y, .x)) %>%
  pivot_longer(-c(iter, parameter), names_to = 'fg') %>%
  pivot_wider(id_cols = c(iter, fg), names_from = parameter, values_from = value) %>%
  mutate(fg = as.numeric(substr(fg, 2, 2)))

# Get fitted values for each iteration (and FG)
fitted_linear_diam <- par_draws_linear_diam %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = dbh_pred, y_fit = threeseg_log(x = dbh_pred, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))

fitted_mixed_diam <- par_draws_mixed_diam %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = dbh_pred, y_fit = threeseg_log(x = dbh_pred, alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)))

# Get quantiles of the fitted values
df_col_names <- c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')

fittedquant_linear_diam <- fitted_linear_diam %>%
  group_by(x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('dbh', df_col_names))

fittedquant_mixed_diam <- fitted_mixed_diam %>%
  group_by(fg, x) %>%
  group_modify(~ data.frame(t(quantile(.$y_fit, probs = qprobs)))) %>%
  setNames(c('fg', 'dbh', df_col_names))

# correction factor for log-log regression plot
# =============================================

# Get fitted values for every data point at every iteration (posterior estimates of linear predictor)
linpred_linear_diam <- par_draws_linear_diam %>%
  mutate(iter = 1:nrow(.)) %>%
  group_by(iter) %>%
  group_modify(~ data.frame(x = dat_all$x, y = dat_all$y, y_fit = threeseg_log(x = dat_all$x, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))

linpred_mixed_diam <- par_draws_mixed_diam %>% 
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = dat_fg$x[dat_fg$fg == .$fg[1]],
                            y = dat_fg$y[dat_fg$fg == .$fg[1]],
                            y_fit = threeseg_log(x =  dat_fg$x[dat_fg$fg == .$fg[1]], alpha = .$coef_alpha, beta = c(.$coef_beta_low, .$coef_beta_mid, .$coef_beta_high), tau_low = .$coef_tau_low, tau_high = .$coef_tau_high)),
               keep = TRUE)

# Take log of fitted values and get residuals by subtracting log of observed y from log of fitted y
linpred_linear_diam <- linpred_linear_diam %>%
  mutate(resid = log(y_fit) - log(y))

linpred_mixed_diam <- linpred_mixed_diam %>%
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

cf_mixed_diam <- linpred_mixed_diam %>%
  group_by(iter) %>%
  summarize(ssq_resid = sum(resid^2),
            n = n()) %>%
  mutate(see = (ssq_resid / (n - k_mixedmodel))^0.5,
         cf = exp((see^2)/2))

# Extract quantiles of correction factor
cfquant_linear_diam <- as.data.frame(t(quantile(cf_linear_diam$cf, probs = qprobs))) %>%
  setNames(df_col_names)

cfquant_mixed_diam <- as.data.frame(t(quantile(cf_mixed_diam$cf, probs = qprobs))) %>%
  setNames(df_col_names)

# Write the extracted output to files
# ===================================

# Write parameters and CI (transform taus to regular units, not log units)
params_3seg_diam <- params_3seg_diam %>%
  mutate(fg = factor(fg, levels = c('all', 1:5), labels = c('all', 'fast', 'large pioneer', 'slow', 'small breeder', 'medium')),
         parameter = factor(parameter, levels = c('alpha', 'beta_low', 'beta_mid', 'beta_high', 'tau_low', 'tau_high'))) %>%
  arrange(fg, parameter) %>%
  mutate_at(vars('mean', starts_with('q')), ~ if_else(grepl('tau', parameter), 10 ^ ., .)) %>%
  mutate(year = 1995, model = '3 segment') %>%
  select(year, model, everything())

write_csv(params_3seg_diam, 'data/data_piecewisefits/richnessbydiameter_paramci_by_fg.csv')

# Write fitted values and CI
fittedquant_diam <- data.frame(
  year = 1995,
  model = '3 segment',
  bind_rows(data.frame(fg = 'all', fittedquant_linear_diam), 
                              fittedquant_mixed_diam %>% ungroup %>% mutate(fg = paste0('fg', fg))))

write_csv(fittedquant_diam, 'data/data_piecewisefits/richnessbydiameter_ci_by_fg.csv')

# Write correction factors
cfquant_diam <- data.frame(fg = c('all', 'by_fg'), bind_rows(cfquant_linear_diam, cfquant_mixed_diam))

write_csv(cfquant_diam, 'data/data_piecewisefits/richnessbydiameter_cf_by_fg.csv')

# Create clean summary tables ---------------------------------------------

# Only need to clean up the parameter tables.
params_3seg_diam %>%
  mutate(parameter = factor(parameter, labels = c('intercept', 'slope small trees', 'slope midsize trees', 'slope large trees', 'cutoff small to midsize', 'cutoff midsize to large'))) %>%
  write_csv('data/clean_summary_tables/clean_parameters_richnessbydiameter.csv')

# Create plotting data ----------------------------------------------------

# Observed data
obs_richnessbydiameter <- data.frame(year = 1995, bind_rows(bin_all, bin_x_fg)) %>%
  mutate(abundance_by_bin_width = n_individuals / (bin_max - bin_min))

# Fitted data (the same as the ci data frame.)

write_csv(obs_richnessbydiameter, 'data/data_forplotting/obs_richnessbydiameter.csv')
write_csv(fittedquant_diam, 'data/data_forplotting/fitted_richnessbydiameter.csv')

# Ratio trends ------------------------------------------------------------

# Modified 27 June 2022: this is now done for all six pairwise combinations of the four named FGs instead of just selected ones.
# For the diameter fits and light fits.
# We divide the fitted values for fast/slow and pioneer/breeder to get the ratios
# The two correction factors cancel out in numerator and denominator, so no need to multiply out by it.

fg_pairs <- combn(1:4, 2, simplify = FALSE)
fg_titles <- c('fast', 'tall', 'slow', 'short')

# Fitted values of ratios (raw)
ratio_diam <- map(fg_pairs, function(fgs) {
  fitted_mixed_diam %>%
    filter(fg %in% fgs) %>%
    pivot_wider(id_cols = c(iter, x), names_from = fg, values_from = y_fit) %>%
    mutate(ratio = .data[[as.character(fgs[1])]]/.data[[as.character(fgs[2])]]) 
})

# Quantiles of fitted values of ratios
ratioquant_diam <- map(ratio_diam, function(dat) {
  dat %>%
    group_by(x) %>%
    group_modify(~ as.data.frame(t(quantile(.$ratio, probs = qprobs)))) %>%
    setNames(c('dbh', df_col_names))
})

# Fitted slopes of ratios (raw)
log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)  

fittedslope_mixed_diam <- fitted_mixed_diam %>%
  group_by(iter, fg) %>%
  group_modify(~ data.frame(x = exp(midpts(log(.$x))), slope = log_slope(.$x, .$y_fit)))

ratioslope_diam <- map(fg_pairs, function(fgs) {
  fittedslope_mixed_diam %>%
    filter(fg %in% fgs) %>%
    pivot_wider(id_cols = c(iter, x), names_from = fg, values_from = slope) %>%
    mutate(slope = .data[[as.character(fgs[1])]] - .data[[as.character(fgs[2])]]) 
})

# Quantiles of fitted slopes of ratios
ratioslopequant_diam <- map(ratioslope_diam, function(dat) {
  dat %>%
    group_by(x) %>%
    group_modify(~ as.data.frame(t(quantile(.$slope, probs = qprobs)))) %>%
    setNames(c('dbh', df_col_names))
})

# Write ratio trends to file
ratio_names <- map(combn(fg_titles, 2, simplify = FALSE), paste, collapse = ':')
ratio_diam_table <- map2_dfr(ratioquant_diam, ratio_names, ~ data.frame(ratio = .y, .x))

write_csv(ratio_diam_table, 'data/data_piecewisefits/ratio_fittedvalues_richnessbydiameter_allpairwise.csv')

# Also put a copy in the data_forplotting directory
write_csv(ratio_diam_table, 'data/data_forplotting/fitted_richnessbydiameter_ratio_allpairwise.csv')

# Observed richness bin values
obs_richness_ratio_diam <- map2_dfr(fg_pairs, ratio_names, function(fgs, ratio_name) {
  bin_top <- bin_x_fg %>% filter(fg == paste0('fg', fgs[1]))
  bin_bottom <- bin_x_fg %>% filter(fg == paste0('fg', fgs[2]))
  tibble(
    ratio_name = ratio_name,
    bin = bin_top$bin, bin_midpoint = bin_top$bin_midpoint, bin_min = bin_top$bin_min, bin_max = bin_top$bin_max,
    richness_ratio = bin_top$richness / bin_bottom$richness,
    lowest_n = pmin(bin_top$n_individuals, bin_bottom$n_individuals),
    lowest_rich = pmin(bin_top$richness, bin_bottom$richness)
  ) %>%
    mutate_if(is.double, ~ if_else(is.finite(.), ., as.numeric(NA)))
})

write_csv(obs_richness_ratio_diam, 'data/data_forplotting/obs_richnessbydiameter_ratio_allpairwise.csv')

# Create clean data frame of fitted slope quantiles, with one slope per segment
# Extract slope at midpoints between cutoffs to get this value

# Get cutoff points
fg_names_old <- c('fast', 'large pioneer', 'slow', 'small breeder')

cutoffs_diam <- map(fg_pairs, function(fgs) {
  params_3seg_diam %>%
    filter(fg %in% fg_names_old[fgs], grepl('tau', parameter)) %>%
    pull(q50) %>% sort
})

# Extract fitted slope values as close to the midpoint as possible
ratioslopesegment_diam <- map2(ratioslopequant_diam, cutoffs_diam, function(slopes, cutoffs) {
  slopes %>%
    ungroup %>%
    mutate(segment = cut(dbh, breaks = c(1, cutoffs[-1], 285), include.lowest = TRUE)) %>%
    group_by(segment) %>%
    filter(abs(dbh-median(dbh)) == min(abs(dbh-median(dbh)))) %>%
    slice(1) %>%
    filter(!is.na(segment)) %>%
    select(-dbh) %>%
    select(segment, everything())
})

ratioslope_all_clean <- 
  map2_dfr(ratioslopesegment_diam, ratio_names, ~ tibble(ratio_name = .y, scaling_variable = 'diameter', .x))


write_csv(ratioslope_all_clean, 'data/clean_summary_tables/clean_richness_ratio_fitted_slopes_allpairwise.csv')


# Abundance vs richness mixed effects model -------------------------------

# For diameter fit.
# Use brms because it is a straightforward mixed effects model that does not require custom Stan code.

library(brms)

bin_x_fg <- bin_x_fg %>% mutate(abundance_by_bin_width = n_individuals / (bin_max - bin_min))

fit_richxabund_diam <- brm(log_richness ~ log_abundance + (log_abundance|fg),
                           data = bin_x_fg %>% 
                             filter(!fg %in% 'unclassified', richness > 0) %>%
                             mutate(log_richness = log10(richness_by_bin_width),
                                    log_abundance = log10(abundance_by_bin_width)),
                           chains = 3, iter = 10000, warmup = 9000, seed = 123)

# Extract coefficients and fitted values
# Use global mean and intercept (fixef) to get fitted line for all trees.

# coefficients (clean up)
coef_richxabund_diam <- coef(fit_richxabund_diam, probs = qprobs)
fixef_richxabund_diam <- fixef(fit_richxabund_diam, probs = qprobs)

params_richxabund_diam <- data.frame(fg = c('all', 'fast', 'large pioneer', 'slow', 'small breeder', 'medium'),
                                     parameter = rep(c('intercept', 'slope'), each = 6),
                                       rbind(fixef_richxabund_diam[1,],
                                             coef_richxabund_diam[[1]][,,1], 
                                             fixef_richxabund_diam[2,],
                                             coef_richxabund_diam[[1]][,,2])) %>%
  setNames(c('fg', 'parameter', 'mean', 'sd', df_col_names))

predict_dat <- expand_grid(fg = c(NA, paste0('fg', 1:5)), 
                           log_abundance = seq(-2, 5, length.out = 500))


fittedraw_richxabund_diam <- fitted(fit_richxabund_diam, newdata = predict_dat, summary = FALSE)
fittedquant_richxabund_diam <- apply(10^fittedraw_richxabund_diam, 2, quantile, probs = qprobs)

fitted_richxabund_diam <- data.frame(abundance_by_bin_width = 10^predict_dat$log_abundance,
                                     fg = predict_dat$fg,
                                     t(fittedquant_richxabund_diam)) %>%
  setNames(c('abundance_by_bin_width', 'fg', df_col_names)) %>%
  mutate(fg = if_else(is.na(fg), 'all', fg))

write_csv(params_richxabund_diam, 'data/clean_summary_tables/clean_parameters_richnessvsabundance.csv')
write_csv(fitted_richxabund_diam, 'data/data_forplotting/fitted_richnessvsabundance.csv')
