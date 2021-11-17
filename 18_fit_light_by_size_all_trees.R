# Script 18
# Get slopes of light per crown area / diameter, light per crown volume / diameter, and light per leaf area / diameter

library(tidyverse)
library(brms)

# Load data
load('ForestLight/data/rawdataobj_alternativecluster.r')

# Do log linear regression of light received / crown area versus dbh

regdata <- alltree_light_95 %>%
  mutate(light_captured_leafarea = light_captured/leaf_area) %>%
  select(dbh_corr, light_received_byarea, light_received_byvolume, light_captured_leafarea) %>%
  setNames(c('dbh', 'light_area', 'light_volume', 'light_captured_leafarea'))

# Plot for sanity check
#ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10()

# fit model.
options(mc.cores = 3)
set.seed(111)

reg_area <- brm(log10(light_area) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)
reg_volume <- brm(log10(light_volume) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)
reg_leafarea <- brm(log10(light_captured_leafarea) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)

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

leafarea_log_fitted <- cbind(1, log10(regdata$dbh)) %*% fixef(reg_leafarea)[,'Estimate'] # Fitted values
leafarea_residuals <- log10(regdata$light_captured_leafarea) - leafarea_log_fitted
leafarea_var_fitted <- var(leafarea_log_fitted)
leafarea_var_resid <- var(leafarea_residuals)
leafarea_rsq <- leafarea_var_fitted / (leafarea_var_fitted + leafarea_var_resid)

# Area 0.668, Volume 0.258, Leaf Area 0.592


# save summary statistics
coef_table <-rbind(
  bind_rows(data.frame(regression = 'light per unit crown area versus dbh', parameter = c('intercept', 'slope'), fixef(reg_area)),
            data.frame(regression = 'light per unit crown area versus dbh', parameter = c('r-squared'), Estimate = area_rsq)),
  bind_rows(data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('intercept', 'slope'), fixef(reg_volume)),
            data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('r-squared'), Estimate = volume_rsq)),
  bind_rows(data.frame(regression = 'light captured per unit leaf area versus dbh', parameter = c('intercept', 'slope'), fixef(reg_leafarea)),
            data.frame(regression = 'light captured per unit leaf area versus dbh', parameter = c('r-squared'), Estimate = leafarea_rsq))
  )

write_csv(coef_table, 'ForestLight/data/clean_summary_tables/clean_parameters_lightperareaandvolumebydiameter.csv')

# save fitted values and credible intervals
dbh_pred <- exp(seq(log(1), log(315), length.out = 101))

fitted_area <- fitted(reg_area, newdata = data.frame(dbh = dbh_pred), summary = TRUE)
fitted_volume <- fitted(reg_volume, newdata = data.frame(dbh = dbh_pred), summary = TRUE)
fitted_leafarea <- fitted(reg_leafarea, newdata = data.frame(dbh = dbh_pred), summary = TRUE)

fitted_all <- data.frame(dbh = dbh_pred,
                         rbind(data.frame(fit = 'light per area', fitted_area),
                               data.frame(fit = 'light per volume', fitted_volume),
                               data.frame(fit = 'light captured per leaf area', fitted_leafarea))) %>%
  rename(q50 = Estimate, q025 = Q2.5, q975 = Q97.5) %>%
  mutate(q50 = 10^q50, q025 = 10^q025, q975 = 10^q975) %>%
  select(fit, dbh, q025, q50, q975)

write_csv(fitted_all, 'ForestLight/data/data_forplotting/fitted_lightbysizealltrees_fig1.csv')
