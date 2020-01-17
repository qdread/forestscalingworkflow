# Allometric coefficients by species
# QDR / Forestlight / 11 Nov 2019

# For tree height and crown area use Martinez-Cano et al. 2019 Table S2 (Michaelis-Menten for height, power function for crown area)
# For crown depth, use raw data provided by Stephanie Bohlman.

library(tidyverse)
library(forestscaling)
library(broom)

bohlman_raw <- read_csv('data/allometry_bci_trees_bohlman.csv')
martinezS2 <- read_csv('/data/allometry_martinez-cano_table_s2.csv')


table(bohlman_raw$species)
length(unique(bohlman_raw$species)) # 80 species

# Cross reference allometry with FG lookup.
fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

fg_lookup <- data.frame(species = fgbci$splower, fg = fgbci$fg5)

bohlman_raw <- bohlman_raw %>%
  select(-X9) %>%
  left_join(fg_lookup)

# Check whether the crown area is just a transformation of crown diameter
diffs <- (pi * (bohlman_raw$crowndiam.m/2)^2) - bohlman_raw$crownarea.m2
range(diffs, na.rm=TRUE) # Yes, they are the same within tolerance.

# Calculate allometries for crown depth from Bohlman's raw data.
fit_and_r2 <- function(dat) {
  mod <- lm(log(value) ~ log(dbh.cm), data = dat)
  data.frame(tidy(mod, conf.int = TRUE), r2 = summary(mod)$r.sq)
}

# Species level coefficients
depth_coefs_sp <- bohlman_raw %>%
  mutate(value = crowndepth.m) %>%
  group_by(fg, species) %>%
  group_modify(~ fit_and_r2(.), keep = TRUE)

# Also calculate group-level coefficients to be used for the species without allometry or the ones with insufficient precision.
depth_coefs_fg <- bohlman_raw %>% 
  mutate(value = crowndepth.m) %>%
  group_by(fg) %>%
  group_modify(~ fit_and_r2(.), keep = TRUE)

# For the species that have neither allo nor fg, need to do one allometry for all.
depth_coefs_overall <- bohlman_raw %>% 
  mutate(value = crowndepth.m) %>%
  group_modify(~ fit_and_r2(.), keep = TRUE)

# Number of individuals per species
allo_n <- bohlman_raw %>% 
  group_by(species) %>%
  summarize(n = sum(!is.na(crowndepth.m)))

depth_coefs_sp <- depth_coefs_sp %>% left_join(allo_n)

# Inspect extreme values of R2 and slope
depth_coefs_sp %>% filter(!term %in% "(Intercept)") %>% arrange(r2) %>% print(n=25)
depth_coefs_sp %>% filter(!term %in% "(Intercept)") %>% arrange(estimate)


# Get coefficients for all species ----------------------------------------

r2_cutoff <- 0.25 # Do not use any with R squared less than 0.25. (will get rid of 7 species)

# Each species needs 8 coefficients, slope and intercept for each of the four allometric equations.

# Make very wide DF with each row a species
depth_coefs_sp_wide <- depth_coefs_sp %>%
  select(-std.error, -statistic, -p.value, -conf.low, -conf.high, -n) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`, slope = `log(dbh.cm)`) %>%
  mutate(slope = if_else(r2 > r2_cutoff, slope, as.double(NA)),
         intercept = if_else(r2 > r2_cutoff, intercept, as.double(NA)))

# Make wide DF with each row a functional group
depth_coefs_fg_wide <- depth_coefs_fg %>%
  select(fg, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`, slope = `log(dbh.cm)`) %>%
  filter(!is.na(fg))

# Also make wide DF for all species together, to use for the unclassified species
depth_coefs_overall_wide <- depth_coefs_overall %>%
  select(term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`, slope = `log(dbh.cm)`) %>%
  mutate(species = 'other')

# Join individual species coefficients with the 282 species lookup table
depth_coefs_allspecies <- fg_lookup %>%
  left_join(depth_coefs_sp_wide)

# Make fg level table and species level table long so they can be joined
depth_coefs_fg_long <- depth_coefs_fg_wide %>%
  pivot_longer(-fg) %>%
  rename(value_fg = value)

depth_coefs_allspecies_long <- depth_coefs_allspecies %>%
  select(-r2) %>%
  pivot_longer(-c(fg, species)) %>%
  left_join(depth_coefs_fg_long)

depth_coefs_allspecies_filled <- depth_coefs_allspecies_long %>%
  mutate(value = if_else(is.na(value), value_fg, value)) %>%
  select(-value_fg) %>%
  pivot_wider(names_from = name, values_from = value)

depth_coefs_allspecies_filled <- bind_rows(depth_coefs_allspecies_filled, depth_coefs_overall_wide) %>%
  rename(crowndepth_a = intercept,
         crowndepth_b = slope)

# Add correction factor for log-transformation bias
mod_all <- lm(log(crowndepth.m) ~ log(dbh.cm), data = bohlman_raw)
depth_corr_factor <- log_correction_factor(log_y = log(bohlman_raw$crowndepth.m), log_y_pred = predict(mod_all), n_pars = 2)

# Add Martinez coefficients -----------------------------------------------

height_area_coefs <- martinezS2 %>%
  select(species, height_a, height_b, height_k, area_a, area_b)

all_coefs <- left_join(depth_coefs_allspecies_filled, height_area_coefs)

# fill in hard-coded coefficients for all species without allometries for height and area
# Use same relationship regardless of functional group, because the raw data are not available to calculate relationship by FG 
# (but all-species allometry is available, included in Table S3 for height, gMM without covariates, and table S4 for area, power fn without covariates)
height_allspecies_params <- c(57.17, 0.7278, 21.57)
area_allspecies_params <- c(0.5659, 1.341)

no_allo_rows <- which(is.na(all_coefs$height_a))
all_coefs[no_allo_rows, c('height_a','height_b','height_k','area_a','area_b')] <- t(replicate(length(no_allo_rows), c(height_allspecies_params, area_allspecies_params)))

# Also fill in hard-coded correction factors for height and for area.
all_coefs$height_corr_factor <- 1.017
all_coefs$area_corr_factor <- 1.163
all_coefs$crowndepth_corr_factor <- depth_corr_factor

# Write to CSV
write_csv(all_coefs, 'data/allometry_final.csv')


