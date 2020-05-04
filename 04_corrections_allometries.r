# Script 4: Correct growth increments, remove outliers, join light data, and apply allometries.

# modified 14 April 2020: use Kitajima et al. 2005 data to establish scaling relationship between LAI and crown depth.
# modified 09 January 2020: add column flagging new recruits so that they can be excluded from later analyses
# modified 11 November 2019: replace overall allometries with species-specific allometries
# modified 30 October 2019: move binning to another script.
# modified 22 October 2019 to correctly add new recruits to production -- also make the code a lot tidier.
# modified 25 June 2018 to use Nadja's newest FG classification

# Each bin dataframe should include the following
# 1 all individuals
# 2 all individuals classified in a FG
# 3-7 FGs 1 through 5
# 8 individuals not classified in a FG

library(tidyverse)
library(forestscaling)

group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')

load('data/bciqcrun.RData')

# Create lists over which to map functions.
bci_full <- mget(paste0('bci.full', 1:7))
bci_stem <- mget(paste0('bci.stem', 1:7))

# Get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

bci_full <- map(bci_full, ~ filter(., !young & !edge))
bci_stem <- map(bci_stem, ~ filter(., !young & !edge))

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

# Load Nadja's data (new functional groups 25 June)
# fg5 is the new column (we originally used fg from the older df)
fgbci <- read.table('data/fgroups_dynamics_new.txt', stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Currently X1new is high for slow species and low for fast species
# Currently X2new is high for pioneer species and low for breeder species
# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new

# Read species-level allometric coefficients
all_coefs <- read_csv('data/allometry_final.csv')

# 2. Convert all census data frames to the correct units (dbh in cm, not mm; agb in kg, not Mg)
# Simultaneously, calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.
# Use harmonic mean to annualize the 5-year biomass increment into a 1-year biomass increment using the small end of the interval
# Also get rid of the young (secondary) forest patches.

# Correction 22 Oct 2019: The old annualizing function returns NaN if the old biomass was 0. 
# if the old biomass was 0, just use 1/5 of the new biomass as the annual growth rate.
# Also annualize to the middle of the interval: 2 to 3
# Also, some of the census intervals are not exactly 5 so we can correct for this with the 

# Get production by taking the difference between successive biomasses and converting to annual rate
bci_production <- map2(bci_full[-7], bci_full[-1], function(old, new) {
  census_interval <- (new$date - old$date)/365.25
  if_else(is.na(old$agb_corr), new$agb_corr / census_interval, annual_increment(meas_old = old$agb_corr, meas_new = new$agb_corr, census_interval = census_interval))
})

# Added 25 Oct 2019: Also get annualized diameter increment (used for model fitting later)
bci_diamgrowthrate <- map2(bci_full[-7], bci_full[-1], function(old, new) {
  census_interval <- (new$date - old$date)/365.25
  if_else(is.na(old$dbh_corr), new$dbh_corr / census_interval, annual_increment(meas_old = old$dbh_corr, meas_new = new$dbh_corr, census_interval = census_interval))
})

# Added 09 Jan 2020: flag the new recruits
bci_recruitstatus <- map2(bci_full[-7], bci_full[-1], function(old, new) {
  is.na(old$dbh_corr) & !is.na(new$dbh_corr)
})

# Amendment 11 Oct. 2017 : Add dbh increment (raw), used to detect production outliers.

bci_dbhincs <- map2(bci_full[-7], bci_full[-1], function(old, new) {
  dbh_old <- if_else(is.na(old$dbh_corr), 0, old$dbh_corr) 
  return(data.frame(dbhinc = new$dbh_corr - dbh_old, dbhlastcensus = dbh_old))
})

# Add columns to the bci_full data frames.

bcicensusdat <- map(1:6, ~ cbind(bci_full[-1][[.]], production = bci_production[[.]], diam_growth_rate = bci_diamgrowthrate[[.]], bci_dbhincs[[.]], recruit = bci_recruitstatus[[.]]))
bcicensusdat <- map(bcicensusdat, function(x) x %>%
                      filter(DFstatus == 'alive') %>%
                      mutate(dbh_corr = dbh_corr/10,
                             agb_corr = agb_corr * 1000,
                             production = production * 1000,
                             diam_growth_rate = diam_growth_rate / 10,
                             dbhinc = dbhinc/10,
                             dbhlastcensus = dbhlastcensus/10) %>%
                      filter(!young) )

####
# Amendment 11 Oct.
# Remove production outliers.
# Criteria: 
# 1. If the tree was not in the previous census and its current dbh is >10, remove.
# 2. If the tree grew over 20 cm dbh between censuses, remove.

bcicensusdat <- map(bcicensusdat, function(x) {
  x$production[x$dbhinc > 20 | (x$dbhinc > 10 & x$dbhlastcensus == 0)] <- NA
  x
})


# 3. Join with shade tolerance (Rüger) groups.
# Here is where the fg5 is assigned as the fg classification we are using.

fgbci_less <- fgbci %>%
  select(sp, grform, fg5, PC_slow_to_fast, PC_breeder_to_pioneer) %>%
  rename(fg = fg5, X1 = PC_slow_to_fast, X2 = PC_breeder_to_pioneer) %>%
  mutate(sp = tolower(sp))

bcicensusdat <- map(bcicensusdat, ~ left_join(., fgbci_less))

# 4. Join the appropriate years with Nadja's light data.
# bcicensusdat[[2]] is 1985-1990
# bcicensusdat[[3]] is 1990-1995

bcicensusdat[[2]] <- bcicensusdat[[2]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth8590 %>% select(tag, light, young, edge))

bcicensusdat[[3]] <- bcicensusdat[[3]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth9095 %>% select(tag, light, young, edge))

# 5. Run allometries to get total light received.

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# Get tree height and crown dimensions from dbh
# Correction made 03 Oct 2019: correct so that we are dividing by volume, not multiplying
# Added 25 Oct 2019: include crown depth so that we can correct for incident light capture percentage using light extinction coefficient.
# Added 28 Oct 2019: include the function to estimate percent light captured given light extinction coefficient of 0.5
# Added 11 Nov 2019: Use species specific parameters generated in previous scripts

overall_k <- 0.5 # Roughly the mean light extinction coefficient for the Panamanian species in Kitajima et al. 2005.

# Get slope and intercept relating LAI to crown depth, from Kitajima et al. 2005 data.
library(brms)

lai_depth <- read_csv('data/LAI_Depth_Kitajima2005.csv') %>%
  mutate(log_LAI = log10(LAI),
         log_depth = log10(Depth)) %>%
  filter(!Species %in% 'Cecropia')

# Fit Bayesian random-intercept, random-slope model
lai_fit <- brm(log_LAI ~ log_depth + (log_depth | Species), data = lai_depth, family = 'gaussian', chains = 3, iter = 15000, warmup = 12500, seed = 111)
# Extract parameter estimates for slope and intercept
lai_coefs <- summary(lai_fit)$fixed[,'Estimate']

for (i in 2:3) {
  
  # Generate lookup table for coefficients
  coef_table <- bcicensusdat[[i]] %>%
    transmute(species = sp) %>%
    left_join(all_coefs)
  # fill in "other"
  coef_table[is.na(coef_table$fg),] <- all_coefs[rep(which(all_coefs$species == 'other'), sum(is.na(coef_table$fg))), ]
  
  bcicensusdat[[i]] <- bcicensusdat[[i]] %>%
    mutate(crownarea = coef_table$area_corr_factor * coef_table$area_a * dbh_corr ^ coef_table$area_b,
           crowndepth = coef_table$crowndepth_corr_factor * exp(coef_table$crowndepth_a + coef_table$crowndepth_b * log(dbh_corr)),
           crownvolume = (2/3) * crownarea * crowndepth, # Half-ellipsoid
           light_received = light * crownarea * insol_bci,
           fraction_light_captured = pct_light_captured(depth = 10^lai_coefs[1] * crowndepth^lai_coefs[2], k = overall_k),
           light_received_byarea = light * insol_bci,
           light_received_byvolume = light * fraction_light_captured * crownarea * insol_bci / crownvolume)

}

# Classification of light into 3 groups.
# This is percent of full irradiance, not the total light received.

for (i in 2:3) {
  
  light_groups <- cut(bcicensusdat[[i]]$light, breaks = 3)
  table(light_groups)
  light_groupcodes <- factor(light_groups, 
                             labels = c('Low light','Intermediate light','High light'))
  
  bcicensusdat[[i]]$light_group <- light_groupcodes
  
}

# Create final raw objects for each tree and split by FG.
# fgdat is nested first by FG, then by year.

alltreedat <- map(bcicensusdat, ~ filter(., !is.na(dbh_corr), production > 0))

fgdat <- map(alltreedat, ~ split(., factor(.$fg, exclude = NULL))) %>% transpose # Transposed to ensure compatibility with old code.

# for 2 and 3, make additional data frames that have only trees with light measurements.
alltree_light_90 <- subset(alltreedat[[2]], !is.na(light))
alltree_light_95 <- subset(alltreedat[[3]], !is.na(light))

light_fg_90 <- lapply(fgdat, function(x) subset(x[[2]], !is.na(light)))
light_fg_95 <- lapply(fgdat, function(x) subset(x[[3]], !is.na(light)))


# Save raw data as an object 
save(alltreedat, fgdat, alltree_light_90, alltree_light_95, light_fg_90, light_fg_95, file = 'data/rawdataobj_alternativecluster.r')

# Save raw data for 1995 only, not separated by fg, so that a small object can be versioned on the repo.
save(alltreedat, alltree_light_95, file = 'data/rawdataobj1995.RData')


