# Script to fit to binned richness values, first for all trees then with a mixed model with coefficients for each FG
# Using both diameter and light/area

# The following procedure is done in this script:
# 1. Load data 
# 2. Load and fit models
# 3. Extract output
#   - parameter values
#   - fitted and predicted values
#   - fitted slopes
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

source('R_functions/model_output_extraction_functions.r')

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

# Only extract from three segment model.
# Extract from the all-tree linear model, and the by-group mixed model, for both diameter and light/area.

# Parameter values

params_3seg_diam <- param_values(fit_3seg_mixed_all, )

# Fitted and predicted values
# Fitted slopes
# Bayesian R-squared and correction factor for log-log regression plot

# Create clean summary tables ---------------------------------------------


# Create plotting data ----------------------------------------------------


