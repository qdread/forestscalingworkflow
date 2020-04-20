# Script 18
# Get slopes of light per crown area / diameter, and light per crown volume / diameter

library(tidyverse)
library(brms)

# Load data
load('data/rawdataobj_alternativecluster.r')

# Do log linear regression of light received / crown area versus dbh

regdata <- alltree_light_95 %>%
  select(dbh_corr, light_received_byarea, light_received_byvolume) %>%
  setNames(c('dbh', 'light_area', 'light_volume'))

# Plot for sanity check
#ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10()

# fit model.
options(mc.cores = 3)
set.seed(111)

reg_area <- brm(log10(light_area) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)
reg_volume <- brm(log10(light_volume) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)

# Use median parameter estimate to get Bayesian R2.
area_log_fitted <- cbind(1, log10(regdata$dbh)) %*% fixef(reg_area)[,'Estimate'] # Fitted values
area_residuals <- log10(regdata$light_area) - area_log_fitted
area_var_fitted <- var(area_log_fitted)
area_var_resid <- var(area_residuals)
area_rsq <- area_var_fitted / (area_var_fitted + area_var_resid)

volume_log_fitted <- cbind(1, log10(regdata$dbh)) %*% fixef(reg_volume)[,'Estimate'] # Fitted values
volume_residuals <- log10(regdata$light_volume) - volume_log_fitted
volume_var_fitted <- var(volume_log_fitted)
volume_var_resid <- var(volume_residuals)
volume_rsq <- volume_var_fitted / (volume_var_fitted + volume_var_resid)

# Area 0.680, Volume 0.262


# save summary statistics
coef_table <-rbind(
  bind_rows(data.frame(regression = 'light per unit crown area versus dbh', parameter = c('intercept', 'slope'), fixef(reg_area)),
            data.frame(regression = 'light per unit crown area versus dbh', parameter = c('r-squared'), Estimate = area_rsq)),
  bind_rows(data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('intercept', 'slope'), fixef(reg_volume)),
            data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('r-squared'), Estimate = volume_rsq)))

write_csv(coef_table, 'data/clean_summary_tables/fig1_light_by_size_parameters.csv')