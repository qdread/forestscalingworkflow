# Process and bin 1990-1995 mortality data 
# QDR / Forestlight / 03 Oct 2019

# Modified 11 Nov 2019 to use the updated allometries.
# Modified 05 Nov 2019 to use our own package.

# Load data ---------------------------------------------------------------

library(tidyverse)
library(forestscaling)

mort <- read.delim('data/mort_final9095.txt', sep = '\t') # read mortality
load('data/rawdataobj_alternativecluster.r') # load other data, which has the FG assignments
all_coefs <- read_csv('data/allometry_final.csv') # read allometric coefficient data in

# Process data ------------------------------------------------------------

# Join the mortality data with the 1995 tree data so that we know we are using the same individuals as for the other analyses we've done.

# Check 
# How many trees are in 1990 and 1995 dataset put together?
length(union(alltreedat[[2]]$tag, alltreedat[[3]]$tag)) # ~217K all together
length(intersect(alltreedat[[2]]$tag, alltreedat[[3]]$tag)) # ~104K in both (alive in both years)
length(setdiff(alltreedat[[2]]$tag, alltreedat[[3]]$tag)) # ~84K in 1990 but not 1995, so they died
length(setdiff(alltreedat[[3]]$tag, alltreedat[[2]]$tag)) # ~29K in 1995 but not in 1990, so they were born

# Create FG lookup table.
fg_lookup <- unique(rbind(alltreedat[[2]], alltreedat[[3]])[,c('sp','fg')])

# Subset the mort dataset based on the tags included in the dataset we've been using, then join it with the FG lookup table
tags_to_use <- union(alltreedat[[2]]$tag, alltreedat[[3]]$tag)

mort_use <- mort %>%
  filter(tag %in% tags_to_use) %>%
  mutate(sp = tolower(sp)) %>%
  left_join(fg_lookup)

### NOTE: we do not do the dbh correction and agb correction that we do for the other data. It can be done later if necessary.

# Calculate crown area and crown volume for the trees in mort_use, with Bohlman's allometric equations.

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# Added 29 Oct 2019: light extinction coefficient for light received by volume
overall_k <- 0.5 # Roughly the mean light extinction coefficient for the Panamanian species in Kitajima et al. 2005.

# Added 11 Nov 2019: create coefficient lookup table
coef_table <- mort_use %>%
  transmute(species = sp) %>%
  left_join(all_coefs)
# fill in "other"
coef_table[is.na(coef_table$fg),] <- all_coefs[rep(which(all_coefs$species == 'other'), sum(is.na(coef_table$fg))), ]

mort_use <- mort_use %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'),
         dbh = dbh/10,
         h = coef_table$height_corr_factor * gMM(x = dbh, a = coef_table$height_a, b = coef_table$height_b, k = coef_table$height_k),    # Height
         crownarea = coef_table$area_corr_factor * coef_table$area_a * dbh ^ coef_table$area_b,
         crowndepth = coef_table$crowndepth_corr_factor * exp(coef_table$crowndepth_a + coef_table$crowndepth_b * log(dbh)),
         crownvolume = (2/3) * crownarea * crowndepth, # Half-ellipsoid
         light_received = light2 * crownarea * insol_bci,
         fraction_light_captured = pct_light_captured(depth = crowndepth, k = overall_k),
         light_received_byarea = light2 * insol_bci,
         light_received_byvolume = light2 * fraction_light_captured * crownarea * insol_bci / crownvolume)

# Do log binning ----------------------------------------------------------

# Create log bins of proportion that died.
# Do this both for light and dbh, and also do this for all trees together and for the individual groups.

### define bin edges to be the same 20 bins used in all other plots.
n_bins <- 20

# Get dbh and light bins from the data.
# Increase maximum value in all cases so that values equal to maximum are counted.
dbh_bin_bounds <- exp(seq(log(min(mort_use$dbh)), log(max(mort_use$dbh) + 0.1), length.out = n_bins + 1))
lightarea_bin_bounds <- exp(seq(log(min(mort_use$light_received_byarea)), log(max(mort_use$light_received_byarea) + 0.1), length.out = n_bins + 1))
lightvolume_bin_bounds <- exp(seq(log(min(mort_use$light_received_byvolume)), log(max(mort_use$light_received_byvolume) + 0.1), length.out = n_bins + 1))

# Get midpoints of bins
dbh_bin_midpoints <- exp(midpts(log(dbh_bin_bounds)))
lightarea_bin_midpoints <- exp(midpts(log(lightarea_bin_bounds)))
lightvolume_bin_midpoints <- exp(midpts(log(lightvolume_bin_bounds)))

### binning by diameter:

# all trees together
bin_mort_dbh_all <- mort_use %>%
  mutate(dbh_bin = as.numeric(cut(dbh, breaks = dbh_bin_bounds, include.lowest = TRUE))) %>%
  group_by(dbh_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
            ) %>%
  ungroup %>%
  mutate(bin_midpoint = dbh_bin_midpoints[dbh_bin]) %>%
  select(bin_midpoint, lived, died, mortality)
  
# by functional group
bin_mort_dbh_fg <- mort_use %>%
  mutate(dbh_bin = as.numeric(cut(dbh, breaks = dbh_bin_bounds, include.lowest = TRUE))) %>%
  group_by(fg, dbh_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
  ) %>%
  ungroup %>%
  mutate(bin_midpoint = dbh_bin_midpoints[dbh_bin]) %>%
  select(fg, bin_midpoint, lived, died, mortality)

bin_mort_dbh_fg <- bind_rows(cbind(fg = 'all', bin_mort_dbh_all), bin_mort_dbh_fg)

### binning by light per area:

# all trees together
bin_mort_lightarea_all <- mort_use %>%
  mutate(lightarea_bin = cut(light_received_byarea, breaks = lightarea_bin_bounds, include.lowest = TRUE)) %>%
  group_by(lightarea_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
  ) %>%
  ungroup %>%
  mutate(bin_midpoint = lightarea_bin_midpoints[lightarea_bin]) %>%
  select(bin_midpoint, lived, died, mortality)

# by functional group
bin_mort_lightarea_fg <- mort_use %>%
  mutate(lightarea_bin = as.numeric(cut(light_received_byarea, breaks = lightarea_bin_bounds, include.lowest = TRUE))) %>%
  group_by(fg, lightarea_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
  ) %>%
  ungroup %>%
  mutate(bin_midpoint = lightarea_bin_midpoints[lightarea_bin]) %>%
  select(fg, bin_midpoint, lived, died, mortality)

bin_mort_lightarea_fg <- bind_rows(cbind(fg = 'all', bin_mort_lightarea_all), bin_mort_lightarea_fg)  

### binning by light per volume:

# all trees together
bin_mort_lightvolume_all <- mort_use %>%
  mutate(lightvolume_bin = cut(light_received_byvolume, breaks = lightvolume_bin_bounds, include.lowest = TRUE)) %>%
  group_by(lightvolume_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
  ) %>%
  ungroup %>%
  mutate(bin_midpoint = lightvolume_bin_midpoints[lightvolume_bin]) %>%
  select(bin_midpoint, lived, died, mortality)

# by functional group
bin_mort_lightvolume_fg <- mort_use %>%
  mutate(lightvolume_bin = as.numeric(cut(light_received_byvolume, breaks = lightvolume_bin_bounds, include.lowest = TRUE))) %>%
  group_by(fg, lightvolume_bin) %>%
  summarize(lived = sum(alive == 1),
            died = sum(alive == 0),
            mortality = died / (lived + died)
  ) %>%
  ungroup %>%
  mutate(bin_midpoint = lightvolume_bin_midpoints[lightvolume_bin]) %>%
  select(fg, bin_midpoint, lived, died, mortality)

bin_mort_lightvolume_fg <- bind_rows(cbind(fg = 'all', bin_mort_lightvolume_all), bin_mort_lightvolume_fg)  


# Write results -----------------------------------------------------------

# Combine into one df

bin_mort <- bind_rows(cbind(variable = 'dbh', bin_mort_dbh_fg),
                      cbind(variable = 'light_per_area', bin_mort_lightarea_fg),
                      cbind(variable = 'light_per_volume', bin_mort_lightvolume_fg))

write_csv(bin_mort, 'data/data_forplotting/obs_mortalitybins.csv')
write_csv(mort_use, 'data/data_forplotting/obs_mortalityindividuals.csv')
