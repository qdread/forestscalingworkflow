# Fit scaling relationships for FG1,2,3,4, 
# individual production vs light/area: log-log regression
# density vs light/area: power law truncated at a higher x-min value.
# These are fits that are used to construct the ratio fits in Figure 6A.
# QDR / 31 March 2020

# Edited 13 April 2020: Include additional computations to create data for plotting, so that the plotting script only contains plotting code.

library(tidyverse)
library(rstan)
library(forestscaling)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

# Load data
load('data/rawdataobj_alternativecluster.r') # doesn't include imputed values

# What is the truncation point
table(round(alltree_light_95$light_received_byarea))[1:20] # Use 7 as the number since it is the mode of the distribution.
lightbin <- logbin(alltree_light_95$light_received_byarea, n = 20)
ggplot(lightbin, aes(x = bin_midpoint, y = bin_value)) + geom_point() + scale_x_log10() + scale_y_log10() # Still looks like 7.

# Create data objects
x_min <- 7

get_stan_data <- function(dat, x_min) with(dat, list(N = nrow(dat), x = dat$light_received_byarea, y = dat$production, x_min = x_min))

stan_data_list <- alltree_light_95 %>%
  filter(!recruit) %>%
  filter(!fg %in% 5, !is.na(fg), light_received_byarea >= x_min) %>%
  mutate(fg = paste0('fg', fg)) %>%
  group_by(fg) %>%
  group_map(~ get_stan_data(., x_min))

# Compile models
mod_dens1 <- stan_model('model_scripts/density1_simplified.stan')
mod_prod1 <- stan_model('model_scripts/production1_nologlik.stan')

# Fit models

# Fitting options
n_chains <- 3
n_iter <- 6000
n_warmup <- 5000

prodfit_alltrees <- map(stan_data_list, ~ sampling(mod_prod1, data = ., chains = n_chains, iter = n_iter, warmup = n_warmup))
densfit_alltrees <- map(stan_data_list, ~ sampling(mod_dens1, data = ., chains = n_chains, iter = n_iter, warmup = n_warmup))

# Pull out the summaries to make sure models look OK
prodfit_summaries <- map(prodfit_alltrees, ~ summary(.)$summary)
densfit_summaries <- map(densfit_alltrees, ~ summary(.)$summary)

# Save fits
save(prodfit_alltrees, densfit_alltrees, prodfit_summaries, densfit_summaries, file = 'data/data_piecewisefits/fits_bylight_forratio.RData')


#### Extract model output to get the fitted values, slopes, etc.
load('data/data_piecewisefits/fits_bylight_forratio.RData')

# source the extra extraction functions that aren't in the package
source('R_functions/model_output_extraction_functions.r')

# Get the statistics on the ratio trends.
la_pred <- logseq(1,412,101)
parnames_prod <- c('beta0', 'beta1')

# Predicted values for each fit.
dens_pred_fg_lightarea <- map(densfit_alltrees, function(fit) {
  pars_fg <- extract(fit, c('alpha')) %>% bind_cols
  pmap(pars_fg, pdf_pareto, x = la_pred, xmin = x_min)
})
prod_pred_fg_lightarea <- map(prodfit_alltrees, function(fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = la_pred)
})

# Get total production including correction factors.

corr_factor <- function(y, y_fit, n_pars) {
  y_fit <- do.call(cbind, y_fit)
  # Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(log(y_fit), 1, log(y))
  
  # Sum of squared residuals
  ssq_resid <- apply(resids^2, 2, sum)
  # Standard error of estimates
  sse <- (ssq_resid / (length(y) - n_pars))^0.5
  # Correction factors
  exp((sse^2)/2)
}

# We need all fitted values for production, not just the 101 values, to calculate correction factor.
prod_pred_all_lightarea <- map2(stan_data_list, prodfit_alltrees, function(dat, fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dat$x)
})

prod_cf_fg_lightarea <- map2(stan_data_list, prod_pred_all_lightarea, ~ corr_factor(y = .x$y, y_fit = .y, n_pars = 2))

totalprod_pred_fg_lightarea <- map2(dens_pred_fg_lightarea, prod_pred_fg_lightarea, ~ do.call(cbind, .x) * do.call(cbind, .y)) %>%
  map2(prod_cf_fg_lightarea, ~ sweep(.x, 2, .y, `*`))

# Take the ratio of the predicted values from each sampling iteration
# Multiply by number of individuals
dens_ratio13_lightarea <- map2(dens_pred_fg_lightarea[[1]], dens_pred_fg_lightarea[[3]], ~ (.x * stan_data_list[[1]]$N) / (.y * stan_data_list[[3]]$N) )
dens_ratio24_lightarea <- map2(dens_pred_fg_lightarea[[2]], dens_pred_fg_lightarea[[4]], ~ (.x * stan_data_list[[2]]$N) / (.y * stan_data_list[[4]]$N))

totalprod_ratio13_lightarea <- totalprod_pred_fg_lightarea[[1]] / totalprod_pred_fg_lightarea[[3]] * (stan_data_list[[1]]$N / stan_data_list[[3]]$N)
totalprod_ratio24_lightarea <- totalprod_pred_fg_lightarea[[2]] / totalprod_pred_fg_lightarea[[4]] * (stan_data_list[[2]]$N / stan_data_list[[4]]$N)

# Generate credible intervals from the ratios
dens_ratio_ci13_lightarea <- do.call(cbind, dens_ratio13_lightarea) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

dens_ratio_ci24_lightarea <- do.call(cbind, dens_ratio24_lightarea) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

totalprod_ratio_ci13_lightarea <- totalprod_ratio13_lightarea %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

totalprod_ratio_ci24_lightarea <- totalprod_ratio24_lightarea %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

# Save fitted values and credible intervals
ratio_fitted_lightarea <- map2_dfr(list(data.frame(ratio = 'fast:slow', variable = 'density'),
                                        data.frame(ratio = 'pioneer:breeder', variable = 'density'),
                                        data.frame(ratio = 'fast:slow', variable = 'total production'),
                                        data.frame(ratio = 'pioneer:breeder', variable = 'total production')),
                                   list(dens_ratio_ci13_lightarea, dens_ratio_ci24_lightarea, totalprod_ratio_ci13_lightarea, totalprod_ratio_ci24_lightarea),
                                   ~data.frame(.x, .y)) %>%
  setNames(c('ratio','variable','q025','q50','q975','light_area'))

write_csv(ratio_fitted_lightarea, 'data/data_piecewisefits/ratio_fittedvalues_lightarea.csv')

# Save parameters in a table
# Parameters for each fit
dens_pars_fg_lightarea <- map(densfit_alltrees, function(fit) {
  extract(fit, c('alpha')) %>% bind_cols
})
prod_pars_fg_lightarea <- map(prodfit_alltrees, function(fit) {
  extract(fit, c('beta0', 'beta1')) %>% bind_cols
})

allpars_lightarea <- map2(dens_pars_fg_lightarea, prod_pars_fg_lightarea, bind_cols) %>%
  map2_dfr(paste0('fg', 1:4), ~ data.frame(fg = .y, .x)) %>%
  mutate(totalprod_slope = beta1 - alpha - 1)

allpars_lightarea_quantiles <- allpars_lightarea %>%
  mutate(id = rep(1:3000, 4)) %>%
  pivot_longer(-c(id, fg)) %>%
  pivot_wider(id_cols = id, names_from = c(fg, name), values_from = value) %>%
  transmute(fast_slow_density_slope = -fg1_alpha + fg3_alpha,
            fast_slow_production_slope = fg1_totalprod_slope - fg3_totalprod_slope,
            pioneer_breeder_density_slope = -fg2_alpha + fg4_alpha,
            pioneer_breeder_production_slope = fg2_totalprod_slope - fg4_totalprod_slope) %>%
  apply(2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) 

data.frame(variable = 'light per area', slope = dimnames(allpars_lightarea_quantiles)[[2]], t(allpars_lightarea_quantiles)) %>%
  setNames(c('variable','ratio','q025', 'q25', 'q50', 'q75', 'q975')) %>%
  write_csv('data/data_piecewisefits/ratio_parameters_lightbyarea.csv')
  
# Code to create plotting data
# ----------------------------

# Get credible intervals --------------------------------------------------

# Multiply density and total production times number of individuals, and divide by area
area_core <- 42.84
fg_names <- c('Fast', 'Tall', 'Slow', 'Short')

dens_pred_fg_lightarea_quantiles <- map2(dens_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  do.call(cbind, dat) %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

dens_pred_dat_lightarea <- map2_dfr(dens_pred_fg_lightarea_quantiles, fg_names,
                          ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

prod_pred_fg_lightarea_quantiles <- map(prod_pred_fg_lightarea, function(dat) {
  do.call(cbind, dat) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
}) 

prod_pred_dat_lightarea <- map2_dfr(prod_pred_fg_lightarea_quantiles, fg_names,
                          ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

totalprod_pred_fg_lightarea_quantiles <- map2(totalprod_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  dat %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

totalprod_pred_dat_lightarea <- map2_dfr(totalprod_pred_fg_lightarea_quantiles, fg_names,
                               ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))


data_to_bin <- alltree_light_95 %>%
  filter(fg %in% 1:4) %>%
  mutate(fg = factor(fg, labels = fg_names)) %>%
  select(fg, light_received_byarea, production)

# Determine bin edges by binning all
binedgedat <- with(data_to_bin, logbin(light_received_byarea, n = 20))

obs_dens_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, edges = binedgedat))

obs_totalprod_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, y = .$production, edges = binedgedat))

obs_indivprod_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ cloudbin_across_years(dat_values = .$production, dat_classes = .$light_received_byarea, edges = binedgedat, n_census = 1))

save(dens_pred_dat_lightarea, prod_pred_dat_lightarea, totalprod_pred_dat_lightarea, obs_dens_lightarea, obs_totalprod_lightarea, obs_indivprod_lightarea, file = 'data/data_forplotting/light_scaling_plotting_data.RData')
