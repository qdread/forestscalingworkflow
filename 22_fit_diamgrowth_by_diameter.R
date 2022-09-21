# Script for additional mortality versus diameter fit
# Hinged piecewise functional form
# Each nonlinear parameter (except hinge smoothness parameter) has functional group as fixed effect
# 21 Sep 2022


# Setup -------------------------------------------------------------------

library(dplyr)
library(readr)
library(brms)
library(tidybayes)
library(ggplot2)
library(tidyr)
library(forestscaling)
library(purrr)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

# Load data
load('data/rawdataobj1995.RData')

diam_data <- alltreedat[[3]] %>%
  filter(!is.na(fg), !recruit) %>%
  select(fg, dbh_corr, diam_growth_rate) %>%
  mutate(fg = paste0('fg', fg))

# Fit model ---------------------------------------------------------------

diam_hinge_fixef_fit <- brm(
  bf(
    log(diam_growth_rate) ~ beta0 + beta1low * (log(dbh_corr) - log(x0)) + (beta1high - beta1low) * delta * log(1 + exp((log(dbh_corr) - log(x0)) / delta)),
    beta0 ~ 0 + fg,
    beta1low ~ 0 + fg,
    beta1high ~ 0 + fg,
    x0 ~ 0 + fg,
    delta ~ 1,
    nl = TRUE
  ),
  data = diam_data, family = gaussian(link = 'identity'),
  prior = c(
    prior(normal(0, 2), nlpar = beta0),
    prior(lognormal(1, 1), nlpar = beta1low, lb = 0),
    prior(lognormal(1, 1), nlpar = beta1high, lb = 0),
    prior(lognormal(1, 1), nlpar = x0, lb = 0),
    prior(exponential(10), nlpar = delta, lb = 0)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 27603,
  file = '~/temp/forestlight/diam_hinge_fixef_brmfit'
)

# Postprocessing and plotting  -------------------------------------

# Set number of bins
numbins <- 20

# Remove trees not belonging to any functional group
alltreedat_classified <- map(alltreedat, ~ filter(., !is.na(fg)))

# Bin classified trees. (log binning of density)
allyeardbh_classified <- map(alltreedat_classified[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)
bin_edges <- c(dbhbin_allclassified$bin_min,dbhbin_allclassified$bin_max[numbins])

qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
diam_bins <- diam_data %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_edges, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  summarize(p = qprobs, q = quantile(diam_growth_rate, probs = qprobs), n = n()) %>%
  pivot_wider(names_from = p, values_from = q) %>%
  setNames(c('fg', 'dbh_bin', 'n', 'q025', 'q25', 'median', 'q75', 'q975')) %>%
  mutate(bin_midpoint = dbhbin_allclassified$bin_midpoint[as.numeric(dbh_bin)])


# Prediction grid: dbh x fg
# Limit ranges to the observed data points and bins with > 20 individuals
diam_max <- diam_bins %>%
  filter(n >= 20) %>%
  group_by(fg) %>%
  filter(bin_midpoint == max(bin_midpoint))

pred_dat <- diam_max %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dbh_corr = exp(seq(log(1), log(.$bin_midpoint), length.out = 101))))

# Get expected values of posterior means (epred) for each combination of diameter and functional group
diamhinge_pred <- pred_dat %>%
  add_epred_draws(diam_hinge_fixef_fit) %>%
  mutate(.epred = exp(.epred))

### Quick diag. plot
ggplot(diam_bins %>% filter(n >= 20), aes(x=bin_midpoint, y = median)) +
  stat_lineribbon(aes(y = .epred, x = dbh_corr), data = diamhinge_pred, size = 0.5) +
  geom_pointrange(aes(ymin = q25, ymax = q75), color = 'gray40') +
  facet_wrap(~ fg) + scale_x_log10(name = 'diameter (cm)') + scale_y_log10(name = 'diameter growth (cm/y)') +
  scale_fill_brewer(palette = 'Blues') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2), strip.background = element_blank())


# Generate parameter tables -----------------------------------------------

diam_hinge_params <- diam_hinge_fixef_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  separate(.variable, into = c('b', 'parameter', 'fg')) %>%
  select(-b) %>%
  mutate(fg = as.integer(substr(fg, 5, 5))) 

diam_hinge_quant <- diam_hinge_params %>%
  group_by(parameter, fg) %>%
  median_qi(.width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(parameter, fg, .value), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, parameter, .value, .lower_0.95, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'parameter', 'median', 'q025', 'q05', 'q25', 'q75', 'q95', 'q975')) %>%
  mutate(parameter = factor(parameter, levels = c('beta0', 'beta1low', 'beta1high', 'x0', 'delta'))) %>%
  arrange(parameter, fg) %>%
  mutate(fg = c(rep(c('fast', 'large pioneer', 'slow', 'small breeder', 'medium'), 4), '(none)'),
         parameter = c(rep(c('beta0 (intercept)', 'beta1low (slope at small size)', 'beta1high (slope at large size)', 'x0 (cutoff between small and large sizes)'), each = 5), 'delta (smoothing parameter for hinge)'))

# Generate plotting data --------------------------------------------------

diam_hinge_plotting_data <- diamhinge_pred %>%
  group_by(fg, dbh_corr) %>%
  median_qi(.epred, .width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(fg, dbh_corr, .epred), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, dbh_corr, .lower_0.95, .lower_0.9, .lower_0.5, .epred, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')) 

write_csv(diam_hinge_quant, 'data/clean_summary_tables/clean_diamgrowthbydiameter_parameters.csv')
write_csv(diam_hinge_plotting_data, 'data/data_forplotting/fitted_diamgrowthbydiameter.csv')
