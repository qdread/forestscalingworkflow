# Light use efficiency

library(tidyverse)
library(forestscaling)

load("ForestLight/data/rawdataobj_withimputedproduction.RData")

# Load light use efficiency extracted data (coefficients and fitted slopes)
load("ForestLight/data/data_piecewisefits/lue_info.RData")


# Binning -----------------------------------------------------------------


alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'),
         LUE = production/light_captured)

# DBH bin edges
numbins <- 20

allyeardbh <- map(alltreedat[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)

# "cloud" bin LUE for all trees and by FG
luebin_all <- with(alltree_light_95, cloudbin_across_years(dat_classes = dbh_corr, dat_values = LUE, edges = dbhbin_all, n_census = 1))
luebin_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  group_modify(~  cloudbin_across_years(dat_classes = .$dbh_corr, dat_values = .$LUE, edges = dbhbin_all, n_census = 1))

luebin <- bind_rows(data.frame(fg = 'all', luebin_all), luebin_fg)


# Turn extracted fitted values into plottable data ------------------------

names(lue_info) <- c('all', paste0('fg', 1:5))
lue_slopes <- map_dfr(lue_info, 'lue_slope_quant', .id = 'fg')
lue_fitted <- map_dfr(lue_info, 'lue_fitted_quant', .id = 'fg')

# Make plot ---------------------------------------------------------------

dat_points <- alltree_light_95 %>% select(fg, dbh_corr, LUE)
dat_points <- bind_rows(dat_points, dat_points %>% mutate(fg = 'all'))

colors <- RColorBrewer::brewer.pal(3, 'RdYlBu')

# Hexbin plot with cloud bins
ggplot() +
  geom_hex(data = dat_points %>% filter(!fg %in% "unclassified"), mapping = aes(x = dbh_corr, y = LUE)) +
  geom_pointrange(data = luebin %>% filter(!fg %in% "unclassified", mean_n_individuals > 20, complete.cases(.)), mapping = aes(x = bin_midpoint, y = median, ymin = q025, ymax = q975)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'diameter (cm)') +
  scale_y_log10(name = parse(text = 'LUE~(kg~y^-1~W^-1)')) +
  scale_fill_gradient2(low = colors[3], mid = colors[2], high = colors[1], trans = 'log10') +
  theme_bw()

# Hexbin plot with fits
ggplot() +
  geom_hex(data = dat_points %>% filter(!fg %in% "unclassified"), mapping = aes(x = dbh_corr, y = LUE)) +
  geom_ribbon(data = lue_fitted, mapping = aes(x = dbh, ymin = q025, ymax = q975), alpha = 0.5, fill = 'black') +
  geom_line(data = lue_fitted, mapping = aes(x = dbh, y = q50)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'diameter (cm)') +
  scale_y_log10(name = parse(text = 'LUE~(kg~y^-1~W^-1)')) +
  scale_fill_gradient2(low = colors[3], mid = colors[2], high = colors[1], trans = 'log10') +
  theme_bw()

# Coefficient plots
ggplot(lue_slopes, aes(x = fg, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_hline(yintercept = 0, linetype = 'dotted', size = 2) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.2) +
  theme_bw() +
  labs(y = 'LUE scaling exponent')
