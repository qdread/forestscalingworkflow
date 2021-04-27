# Script to clean up all new piecewise output and put in tidy CSVs
# QDR / Forestlight / 14 June 2019

# Modified 26 Feb 2020: add summary table for fitted total production slope across intervals
# Modified 20 June: add light params and R2s.
# Modified 10 Sep: include sigma in production fits.

fp_out <- 'clean_summary_tables'

library(tidyverse)

params <- read.csv('finalcsvs/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

density_param_df <- params %>% filter(variable == 'density')
production_param_df <- params %>% filter(variable == 'production')

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')
density_params <- c('alpha', 'alpha_low', 'alpha_high', 'tau', 'alpha_mid', 'tau_low', 'tau_high')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter', 'sigma')
density_params_names <- c('slope', 'slope small trees', 'slope large trees', 'cutoff small to large', 'slope middle trees', 'cutoff small to middle', 'cutoff middle to large')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# Convert parameters in density data frame to slopes.
density_param_df <- density_param_df %>%
  mutate_at(vars(matches('^mean|^q')), ~ case_when(parameter %in% c('alpha', 'alpha_mid', 'alpha_high') ~ -(.x+1),
                                         parameter %in% c('alpha_low') ~ .x-1,
                                         TRUE ~ .x))

# correct to alternative names
production_param_df <- production_param_df %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

density_param_df <- density_param_df %>%
  mutate(parameter_description = density_params_names[match(parameter, density_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment',
                           model == 3 ~ 'three segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

write.csv(production_param_df, file.path(fp_out, 'clean_parameters_production.csv'), row.names = FALSE)
write.csv(density_param_df, file.path(fp_out, 'clean_parameters_density.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                    prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_production.csv'), row.names = FALSE)

# Clean output of informatiion criteria
ics <- read.csv('finalcsvs/piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment',
                      dens_model == 1 ~ 'one segment',
                      dens_model == 2 ~ 'two segment',
                      dens_model == 3 ~ 'three segment'),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup.csv'), row.names = FALSE)


# Individual light output -------------------------------------------------

library(tidyverse)

params <- read.csv('finalcsvs/light_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivlight_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivlight_param_df, file.path(fp_out, 'clean_parameters_individuallight.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/light_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_individuallight.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('finalcsvs/light_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
                      ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_individuallight.csv'), row.names = FALSE)

# Volume output -------------------------------------------------

library(tidyverse)

params <- read.csv('finalcsvs/volume_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
volume_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(volume_param_df, file.path(fp_out, 'clean_parameters_crownvolume.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/volume_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_crownvolume.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('finalcsvs/volume_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_crownvolume.csv'), row.names = FALSE)


# Captured light output -------------------------------------------------

library(tidyverse)

params <- read.csv('finalcsvs/lightcaptured_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivlightcaptured_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivcapturedlight_param_df, file.path(fp_out, 'clean_parameters_individuallightcaptured.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/lightcaptured_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_individuallightcaptured.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('finalcsvs/lightcaptured_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_individuallightcaptured.csv'), row.names = FALSE)

# Leaf area output -------------------------------------------------

library(tidyverse)

params <- read.csv('finalcsvs/leafarea_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
leafarea_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(leafarea_param_df, file.path(fp_out, 'clean_parameters_leafarea.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/leafarea_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_leafarea.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('finalcsvs/leafarea_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_leafarea.csv'), row.names = FALSE)


# Diameter growth output --------------------------------------------------

library(tidyverse)

params <- read.csv('finalcsvs/diamgrowth_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter', 'sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivdiamgrowth_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivdiamgrowth_param_df, file.path(fp_out, 'clean_parameters_individualdiametergrowth.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('finalcsvs/diamgrowth_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_individualdiametergrowth.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('finalcsvs/diamgrowth_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_individualdiametergrowth.csv'), row.names = FALSE)


# Mortality regression ----------------------------------------------------

mort_pars <- read_csv('finalcsvs/mortality_paramci_by_fg.csv')

mort_pars <- filter(mort_pars, !parameter %in% 'lp__', !grepl('fg', parameter))

mort_pars$parameter <- c('global mean intercept', 'intercept std. dev', 'global mean slope', 'slope std.dev', rep('intercept', 5), rep('slope', 5))

mort_pars$fg <- c(rep('--', 4), rep(fg_full_names[1:5], 2))

mort_pars <- mort_pars %>% select(fg, parameter, everything())

names(mort_pars)[6:10] <- c('q025', 'q25', 'q50', 'q75', 'q975')

write_csv(mort_pars, file.path(fp_out, 'clean_parameters_mortality.csv'))


# Production per area by light per area -----------------------------------

# parameters
la_pars <- read_csv('finalcsvs/lightbyarea_paramci_by_fg.csv')

la_par_description <- c('A','b','k','light per area at maximum slope of curve', 'growth per area at maximum slope of curve', 'maximum slope of curve in log-log space')

la_pars_df <- la_pars %>%
  mutate(parameter_description = rep(la_par_description,7),
         fg = fg_full_names[match(fg, fgs)]) %>%
  select(-parameter) %>%
  select(year, fg, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

write_csv(la_pars_df, file.path(fp_out, 'clean_parameters_growthperareabylightperarea.csv'))

# r squared values
r2df <- read.csv('finalcsvs/lightbyarea_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(fg, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_growthperareabylightperarea.csv'), row.names = FALSE)


# Fitted slopes of total production on intervals --------------------------

# Read fitted slopes
fitted_slopes <- read.csv('data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'total_production', prod_model == 1)

# Read parameters to get the cutoff points, reshape to put all in one row
cutoffs <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'density', grepl('tau', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

# Using cutoffs, calculate points at which to extract the fitted slopes
# For model 1 it can be any point, for model 2 midpoint of upper and lower segment, for model 3 midpoint of each of the 3 segments
dbh_range <- c(1, 285)
slope_points <- cutoffs %>%
  group_by(fg) %>%
  mutate(point_1 = exp(mean(log(dbh_range))),
         pointlow_2 = exp(mean(log(c(dbh_range[1], tau)))),
         pointhigh_2 = exp(mean(log(c(dbh_range[2], tau)))),
         pointlow_3 = exp(mean(log(c(dbh_range[1], tau_low)))),
         pointmid_3 = exp(mean(log(c(tau_low, tau_high)))),
         pointhigh_3 = exp(mean(log(c(dbh_range[2], tau_high))))) %>%
  select(-(tau:tau_high)) %>%
  pivot_longer(-fg) %>%
  separate(name, into = c('point', 'dens_model'), sep = '_') %>%
  mutate(dens_model = as.integer(dens_model))

# Join the fitted slope data frame with the cutoff points and extract the slopes at the points closest to the midpoints
slopes_at_points <- fitted_slopes %>%
  left_join(slope_points) %>%
  mutate(diff = abs(dbh - value)) %>%
  group_by(fg, dens_model, point) %>%
  filter(diff == min(diff)) %>%
  select(year, fg, dens_model, point, q025:q975) %>%
  ungroup

# Add median cutoffs to the table
cutoff_medians <- cutoffs %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dens_model = c(1,2,2,3,3,3), 
                            point = c('point','pointlow','pointhigh','pointlow','pointmid','pointhigh'),
                            segment = c('all sized trees', 'small trees', 'large trees', 'small trees', 'midsize trees', 'large trees'),
                            dbh_start = c(1, 1, .$tau, 1, .$tau_low, .$tau_high),
                            dbh_end = c(285, .$tau, 285, .$tau_low, .$tau_high, 285)))

# add descriptive names to table
slope_table_final <- slopes_at_points %>%
  left_join(cutoff_medians %>% ungroup) %>%
  mutate(fg = fg_full_names[match(fg, fgs)],
         dens_model = rep(c('one segment', 'two segment', 'two segment', 'three segment', 'three segment', 'three segment'), times = 7)) %>%
  select(year, fg, dens_model, segment, dbh_start, dbh_end, q025:q975)
  
write.csv(slope_table_final, file.path(fp_out, 'clean_total_production_fitted_slopes.csv'), row.names = FALSE)


# Fitted slopes of total light on intervals -------------------------------

# Read fitted slopes
fitted_slopes <- read.csv('data/data_piecewisefits/light_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'total_incoming_light', prod_model == 1)

# Read parameters to get the cutoff points, reshape to put all in one row
cutoffs <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'density', grepl('tau', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

# Using cutoffs, calculate points at which to extract the fitted slopes
# For model 1 it can be any point, for model 2 midpoint of upper and lower segment, for model 3 midpoint of each of the 3 segments
dbh_range <- c(1, 285)
slope_points <- cutoffs %>%
  group_by(fg) %>%
  mutate(point_1 = exp(mean(log(dbh_range))),
         pointlow_2 = exp(mean(log(c(dbh_range[1], tau)))),
         pointhigh_2 = exp(mean(log(c(dbh_range[2], tau)))),
         pointlow_3 = exp(mean(log(c(dbh_range[1], tau_low)))),
         pointmid_3 = exp(mean(log(c(tau_low, tau_high)))),
         pointhigh_3 = exp(mean(log(c(dbh_range[2], tau_high))))) %>%
  select(-(tau:tau_high)) %>%
  pivot_longer(-fg) %>%
  separate(name, into = c('point', 'dens_model'), sep = '_') %>%
  mutate(dens_model = as.integer(dens_model))

# Join the fitted slope data frame with the cutoff points and extract the slopes at the points closest to the midpoints
slopes_at_points <- fitted_slopes %>%
  left_join(slope_points) %>%
  mutate(diff = abs(dbh - value)) %>%
  group_by(fg, dens_model, point) %>%
  filter(diff == min(diff)) %>%
  select(year, fg, dens_model, point, q025:q975) %>%
  ungroup

# Add median cutoffs to the table
cutoff_medians <- cutoffs %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dens_model = c(1,2,2,3,3,3), 
                            point = c('point','pointlow','pointhigh','pointlow','pointmid','pointhigh'),
                            segment = c('all sized trees', 'small trees', 'large trees', 'small trees', 'midsize trees', 'large trees'),
                            dbh_start = c(1, 1, .$tau, 1, .$tau_low, .$tau_high),
                            dbh_end = c(285, .$tau, 285, .$tau_low, .$tau_high, 285)))

# add descriptive names to table
slope_table_final <- slopes_at_points %>%
  left_join(cutoff_medians %>% ungroup) %>%
  mutate(fg = fg_full_names[match(fg, fgs)],
         dens_model = rep(c('one segment', 'two segment', 'two segment', 'three segment', 'three segment', 'three segment'), times = 7)) %>%
  select(year, fg, dens_model, segment, dbh_start, dbh_end, q025:q975)

write.csv(slope_table_final, file.path(fp_out, 'clean_total_light_fitted_slopes.csv'), row.names = FALSE)

# Fitted slopes of total crown volume on intervals ------------------------

# Read fitted slopes
fitted_slopes <- read.csv('data/data_piecewisefits/volume_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'total_crown_volume', prod_model == 1)

# Read parameters to get the cutoff points, reshape to put all in one row
cutoffs <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'density', grepl('tau', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

# Using cutoffs, calculate points at which to extract the fitted slopes
# For model 1 it can be any point, for model 2 midpoint of upper and lower segment, for model 3 midpoint of each of the 3 segments
dbh_range <- c(1, 285)
slope_points <- cutoffs %>%
  group_by(fg) %>%
  mutate(point_1 = exp(mean(log(dbh_range))),
         pointlow_2 = exp(mean(log(c(dbh_range[1], tau)))),
         pointhigh_2 = exp(mean(log(c(dbh_range[2], tau)))),
         pointlow_3 = exp(mean(log(c(dbh_range[1], tau_low)))),
         pointmid_3 = exp(mean(log(c(tau_low, tau_high)))),
         pointhigh_3 = exp(mean(log(c(dbh_range[2], tau_high))))) %>%
  select(-(tau:tau_high)) %>%
  pivot_longer(-fg) %>%
  separate(name, into = c('point', 'dens_model'), sep = '_') %>%
  mutate(dens_model = as.integer(dens_model))

# Join the fitted slope data frame with the cutoff points and extract the slopes at the points closest to the midpoints
slopes_at_points <- fitted_slopes %>%
  left_join(slope_points) %>%
  mutate(diff = abs(dbh - value)) %>%
  group_by(fg, dens_model, point) %>%
  filter(diff == min(diff)) %>%
  select(year, fg, dens_model, point, q025:q975) %>%
  ungroup

# Add median cutoffs to the table
cutoff_medians <- cutoffs %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dens_model = c(1,2,2,3,3,3), 
                            point = c('point','pointlow','pointhigh','pointlow','pointmid','pointhigh'),
                            segment = c('all sized trees', 'small trees', 'large trees', 'small trees', 'midsize trees', 'large trees'),
                            dbh_start = c(1, 1, .$tau, 1, .$tau_low, .$tau_high),
                            dbh_end = c(285, .$tau, 285, .$tau_low, .$tau_high, 285)))

# add descriptive names to table
slope_table_final <- slopes_at_points %>%
  left_join(cutoff_medians %>% ungroup) %>%
  mutate(fg = fg_full_names[match(fg, fgs)],
         dens_model = rep(c('one segment', 'two segment', 'two segment', 'three segment', 'three segment', 'three segment'), times = 7)) %>%
  select(year, fg, dens_model, segment, dbh_start, dbh_end, q025:q975)

write.csv(slope_table_final, file.path(fp_out, 'clean_total_volume_fitted_slopes.csv'), row.names = FALSE)


# Fitted slopes of total captured light on intervals -------------------------------

# Read fitted slopes
fitted_slopes <- read.csv('data/data_piecewisefits/lightcaptured_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'total_captured_light', prod_model == 1)

# Read parameters to get the cutoff points, reshape to put all in one row
cutoffs <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'density', grepl('tau', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

# Using cutoffs, calculate points at which to extract the fitted slopes
# For model 1 it can be any point, for model 2 midpoint of upper and lower segment, for model 3 midpoint of each of the 3 segments
dbh_range <- c(1, 285)
slope_points <- cutoffs %>%
  group_by(fg) %>%
  mutate(point_1 = exp(mean(log(dbh_range))),
         pointlow_2 = exp(mean(log(c(dbh_range[1], tau)))),
         pointhigh_2 = exp(mean(log(c(dbh_range[2], tau)))),
         pointlow_3 = exp(mean(log(c(dbh_range[1], tau_low)))),
         pointmid_3 = exp(mean(log(c(tau_low, tau_high)))),
         pointhigh_3 = exp(mean(log(c(dbh_range[2], tau_high))))) %>%
  select(-(tau:tau_high)) %>%
  pivot_longer(-fg) %>%
  separate(name, into = c('point', 'dens_model'), sep = '_') %>%
  mutate(dens_model = as.integer(dens_model))

# Join the fitted slope data frame with the cutoff points and extract the slopes at the points closest to the midpoints
slopes_at_points <- fitted_slopes %>%
  left_join(slope_points) %>%
  mutate(diff = abs(dbh - value)) %>%
  group_by(fg, dens_model, point) %>%
  filter(diff == min(diff)) %>%
  select(year, fg, dens_model, point, q025:q975) %>%
  ungroup

# Add median cutoffs to the table
cutoff_medians <- cutoffs %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dens_model = c(1,2,2,3,3,3), 
                            point = c('point','pointlow','pointhigh','pointlow','pointmid','pointhigh'),
                            segment = c('all sized trees', 'small trees', 'large trees', 'small trees', 'midsize trees', 'large trees'),
                            dbh_start = c(1, 1, .$tau, 1, .$tau_low, .$tau_high),
                            dbh_end = c(285, .$tau, 285, .$tau_low, .$tau_high, 285)))

# add descriptive names to table
slope_table_final <- slopes_at_points %>%
  left_join(cutoff_medians %>% ungroup) %>%
  mutate(fg = fg_full_names[match(fg, fgs)],
         dens_model = rep(c('one segment', 'two segment', 'two segment', 'three segment', 'three segment', 'three segment'), times = 7)) %>%
  select(year, fg, dens_model, segment, dbh_start, dbh_end, q025:q975)

write.csv(slope_table_final, file.path(fp_out, 'clean_total_light_captured_fitted_slopes.csv'), row.names = FALSE)

# Fitted slopes of total leaf area on intervals ------------------------

# Read fitted slopes
fitted_slopes <- read.csv('data/data_piecewisefits/leafarea_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'total_leaf_area', prod_model == 1)

# Read parameters to get the cutoff points, reshape to put all in one row
cutoffs <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE) %>%
  filter(variable %in% 'density', grepl('tau', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

# Using cutoffs, calculate points at which to extract the fitted slopes
# For model 1 it can be any point, for model 2 midpoint of upper and lower segment, for model 3 midpoint of each of the 3 segments
dbh_range <- c(1, 285)
slope_points <- cutoffs %>%
  group_by(fg) %>%
  mutate(point_1 = exp(mean(log(dbh_range))),
         pointlow_2 = exp(mean(log(c(dbh_range[1], tau)))),
         pointhigh_2 = exp(mean(log(c(dbh_range[2], tau)))),
         pointlow_3 = exp(mean(log(c(dbh_range[1], tau_low)))),
         pointmid_3 = exp(mean(log(c(tau_low, tau_high)))),
         pointhigh_3 = exp(mean(log(c(dbh_range[2], tau_high))))) %>%
  select(-(tau:tau_high)) %>%
  pivot_longer(-fg) %>%
  separate(name, into = c('point', 'dens_model'), sep = '_') %>%
  mutate(dens_model = as.integer(dens_model))

# Join the fitted slope data frame with the cutoff points and extract the slopes at the points closest to the midpoints
slopes_at_points <- fitted_slopes %>%
  left_join(slope_points) %>%
  mutate(diff = abs(dbh - value)) %>%
  group_by(fg, dens_model, point) %>%
  filter(diff == min(diff)) %>%
  select(year, fg, dens_model, point, q025:q975) %>%
  ungroup

# Add median cutoffs to the table
cutoff_medians <- cutoffs %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dens_model = c(1,2,2,3,3,3), 
                            point = c('point','pointlow','pointhigh','pointlow','pointmid','pointhigh'),
                            segment = c('all sized trees', 'small trees', 'large trees', 'small trees', 'midsize trees', 'large trees'),
                            dbh_start = c(1, 1, .$tau, 1, .$tau_low, .$tau_high),
                            dbh_end = c(285, .$tau, 285, .$tau_low, .$tau_high, 285)))

# add descriptive names to table
slope_table_final <- slopes_at_points %>%
  left_join(cutoff_medians %>% ungroup) %>%
  mutate(fg = fg_full_names[match(fg, fgs)],
         dens_model = rep(c('one segment', 'two segment', 'two segment', 'three segment', 'three segment', 'three segment'), times = 7)) %>%
  select(year, fg, dens_model, segment, dbh_start, dbh_end, q025:q975)

write.csv(slope_table_final, file.path(fp_out, 'clean_total_leaf_area_fitted_slopes.csv'), row.names = FALSE)


# Ratio slope tables ------------------------------------------------------

# Read the data for ratio fitted slopes for diameter and remove duplicated rows
ratioslopes_diameter <- read_csv('data_piecewisefits/ratio_slope_ci.csv')

# Median cutoff points to define segments
cutoffs_fastslow <- c(3.2, 9.8, 42.3, 48.2)
cutoffs_breederpioneer <- c(2.1, 6.0, 18.9, 71.2)


ratioslopes_diameter_fastslow <- ratioslopes_diameter %>%
  filter(ratio == 'fast:slow') %>%
  mutate(segment = cut(dbh, breaks = c(1, cutoffs_fastslow, 285), include.lowest = TRUE)) %>%
  group_by(variable, segment) %>%
  filter(abs(dbh-median(dbh)) == min(abs(dbh-median(dbh)))) %>%
  slice(1) %>%
  filter(!is.na(segment)) %>%
  select(-dbh) %>%
  mutate(scaling_variable = 'diameter') %>%
  select(ratio, scaling_variable, segment, everything())

ratioslopes_diameter_pioneerbreeder <- ratioslopes_diameter %>%
  filter(ratio == 'pioneer:breeder') %>%
  mutate(segment = cut(dbh, breaks = c(1, cutoffs_breederpioneer, 285), include.lowest = TRUE)) %>%
  group_by(variable, segment) %>%
  filter(abs(dbh-median(dbh)) == min(abs(dbh-median(dbh)))) %>%
  slice(1) %>%
  filter(!is.na(segment)) %>%
  select(-dbh) %>%
  mutate(scaling_variable = 'diameter') %>%
  select(ratio, scaling_variable, segment, everything())

ratioslopequantiles_light <- read_csv('data_piecewisefits/ratio_parameters_lightbyarea.csv') %>%
  rename(scaling_variable = variable) %>%
  mutate(segment = 'all', ratio = rep(c('fast:slow','pioneer:breeder'), each = 2), variable = rep(c('density', 'total production'), times = 2))

allratioslopequantiles <- bind_rows(ratioslopes_diameter_fastslow, ratioslopes_diameter_pioneerbreeder, ratioslopequantiles_light) %>%
  select(-q05, -q95) 
  
write.csv(allratioslopequantiles, file.path(fp_out, 'clean_ratio_fitted_slopes.csv'), row.names = FALSE)
