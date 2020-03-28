# Quality Control run on raw BCI data
# Follows Meakem et al. 2017
# QDR 11 July 2017

# Modified 11 Nov. 2019: Replace overall height allometry with species specific generalized Michaelis Menten height allometry.
# Modified 28 Oct. 2019: Also include light data from recruits (not just trees that have a growth rate in 1995)

# Load Condit's BCI data and Nadja's light data.

growth8590 <- read.delim('data/growth_final8590.txt', stringsAsFactors = FALSE)
growth9095 <- read.delim('data/growth_final9095.txt', stringsAsFactors = FALSE)

load('data/recruits1995_light.rdata')

# Combine 90-95 growth+light data with 90-95 recruits' light data

library(MASS)
library(nlme)
library(tidyverse)
library(forestscaling)

recs9095l <- recs9095l %>% mutate(
  tag = as.integer(tag),
  dinc = NA,
  interval = NA,
  sp = toupper(sp)
)

growth9095 <- rbind(growth9095, recs9095l[, names(growth9095)])

# Load census data. 1 = 1980, 2 = 1985, 3 = 1990, 4 = 1995, 5 = 2000, 6 = 2005, 7 = 2010.

for (i in 1:7) {
  load(paste0('data/BCI_raw/bci.full', i, '.rdata'))
  load(paste0('data/BCI_raw/bci.stem', i, '.rdata'))
}

# Get rid of tree ferns and strangler figs
# Tree ferns are already removed
species_to_remove <- c('ficubu','ficuc1','ficuc2','ficuci','ficupe','ficupo')
remove_spp <- function(dat) subset(dat, !sp %in% c(species_to_remove, toupper(species_to_remove)))

for (i in 1:7) {
  assign(paste0('bci.full',i), remove_spp(get(paste0('bci.full',i))))
  assign(paste0('bci.stem',i), remove_spp(get(paste0('bci.stem',i))))
}

growth8590 <- remove_spp(growth8590)
growth9095 <- remove_spp(growth9095)

# Apply taper correction for trees not measured at dbh because of their buttress, returning an estimate of what the true dbh would be.
# See Cushman reference for the formula
# Additional data files helpfully provided by KC Cushman, 13 July 2017.
# Details on WSG measurement protocol are linked on this page: http://www.forestgeo.si.edu/group/Plant+Functional+Traits/

wsg.data <- read.csv('data/wsg.data.csv', stringsAsFactors = FALSE)
taper.parameters <- read.csv('data/taper.parameters.csv', stringsAsFactors = FALSE)
all_coefs <- read.csv('data/allometry_final.csv', stringsAsFactors = FALSE)

# These are the species that need to be corrected for taper.
taper_spp <- c("alsebl", "anacex", "brosal", "cavapl", "ceibpe", "diptpa", "huracr", "other", "quaras", "tab2ar", "tet2pa")        

#    Make a vector of palm species with individuals measured above 1.3 m height (not 
#    included for taper correction)
palmsp <- c("elaeol","sch1zo","socrex")


# Values of the b1 parameter for each of the species that need to be corrected for taper. See Table 3, Cushman et al. 2014.

model1.pars <- taper.parameters[taper.parameters$eqn==1 & taper.parameters$b1 >0,]

#  The best model includes 2010 diameter, 2010 height of measurement, and species group
varmodel <- gls(log(b1)~log(diam.2010)+log(hom.2010)+spb, data=model1.pars, method="ML")
#    Extract coefficients and covariance matrix from the model
varmodel.coefficients <- varmodel$ coefficients
varmodel.cormatrix <- summary(varmodel)$ varBeta
#    Extract coefficients for relationship with 2010 diameter and height of measurement
mdbh <- varmodel$coefficients[2]
mhom <- varmodel$coefficients[3]

# To find b1 for each of the species, we need to calculate it from the equation:
# b1 <- exp(mdbh * log(dbh) + mhom * log(hom) + intercept + species coefficient)

spcoeffs <- rep(0, length(taper_spp))

for (i in 1:length(spcoeffs)) {
  
  # Assign each tree a taper parameter based on its species group, 2010 diameter, and 2010 
  # height of measurement 
     spcoeffs[i] <- ifelse(i == 1,  varmodel.coefficients[1], varmodel.coefficients[1] + varmodel.coefficients[i+2])

}

names(spcoeffs) <- taper_spp



# Apply dbh correction

for (i in 1:7) {
  stem_i <- get(paste0('bci.stem', i))
  full_i <- get(paste0('bci.full', i))
  dbh_corr_stem_i <- taper_correct_dbh(stem_i)
  dbh_corr_full_i <- taper_correct_dbh(full_i)
  assign(paste0('bci.stem', i), cbind(stem_i, dbh_corr_stem_i))
  assign(paste0('bci.full', i), cbind(full_i, dbh_corr_full_i))
}

# Use dbh correction and wsg table to get corrected aboveground biomass values.
# Use species specific height allometries

# Create allometry lookup table for stem and full data
stem_height_coefs <- bci.stem1 %>%
  transmute(species = sp) %>%
  left_join(all_coefs %>% select(species, height_a, height_b, height_k)) %>%
  mutate(height_a = if_else(is.na(height_a), all_coefs$height_a[all_coefs$species %in% 'other'], height_a),
         height_b = if_else(is.na(height_b), all_coefs$height_b[all_coefs$species %in% 'other'], height_b),
         height_k = if_else(is.na(height_k), all_coefs$height_k[all_coefs$species %in% 'other'], height_k))
full_height_coefs <- bci.full1 %>%
  transmute(species = sp) %>%
  left_join(all_coefs %>% select(species, height_a, height_b, height_k)) %>%
  mutate(height_a = if_else(is.na(height_a), all_coefs$height_a[all_coefs$species %in% 'other'], height_a),
         height_b = if_else(is.na(height_b), all_coefs$height_b[all_coefs$species %in% 'other'], height_b),
         height_k = if_else(is.na(height_k), all_coefs$height_k[all_coefs$species %in% 'other'], height_k))


for (i in 1:7) {
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% 
           left_join(wsg.data) %>%
           mutate(height_corr = all_coefs$height_corr_factor[1] * gMM(x = dbh_corr/10, a = stem_height_coefs$height_a, b = stem_height_coefs$height_b, k = stem_height_coefs$height_k),
                  agb_corr = agb.allometry(wsg, dbh_corr/10, height_corr)/1000))
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% 
           left_join(wsg.data) %>%
           mutate(height_corr = all_coefs$height_corr_factor[1] * gMM(x = dbh_corr/10, a = full_height_coefs$height_a, b = full_height_coefs$height_b, k = full_height_coefs$height_k),
                  agb_corr = agb.allometry(wsg, dbh_corr/10, height_corr)/1000))
}



# Remove the patch of secondary forest in the BCI plot
habitat <- read.table('data/habitats2.txt', header = TRUE, stringsAsFactors = FALSE)

# Define area of young forest
young_forest <- subset(habitat, Habitat == 'young')
with(young_forest, plot(GX20, GY20))

# To remove young forest, we need to round all BCI trees to the nearest multiple of 20 lower than their coordinate, and see if it is in one of those.
young_forest_list <- data.frame(t(young_forest[,c('GX20','GY20')]))

# Vectorized function to assign individual stems to young forest area
is_young <- function(x, y) {
  x_quadrat <- plyr::round_any(x, 20, f = floor)
  y_quadrat <- plyr::round_any(y, 20, f = floor)
  quadrat_list <- data.frame(rbind(x_quadrat, y_quadrat))
  !is.na(match(quadrat_list, young_forest_list))
}

for (i in 1:7) {
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% mutate(young = is_young(gx, gy)))
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% mutate(young = is_young(gx, gy)))
}

growth8590 <- growth8590 %>% mutate(young = is_young(gx, gy))
growth9095 <- growth9095 %>% mutate(young = is_young(gx, gy))

# Calculate the area that Nadja has light values for, so we can determine the right denominator for the per-hectare values.
# 20-meter strip around the outside: gx < 20 | gx > 980 | gy < 20 | gy > 480

is_edge <- function(gx, gy) gx < 20 | gx > 980 | gy < 20 | gy > 480

for (i in 1:7) {
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% mutate(edge = is_edge(gx, gy)))
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% mutate(edge = is_edge(gx, gy)))
}

growth8590 <- mutate(growth8590, edge = is_edge(gx, gy))
growth9095 <- mutate(growth9095, edge = is_edge(gx, gy))

# Calculate number of hectares that is neither edge nor young.
# number of square meters in edge
edge_area <- 2 * 20 * 1000 + 2 * 20 * 460
# number of square meters in young forest
young_area <- nrow(young_forest) * 400
# number of square meters in young forest that isn't edge
young_area_core <- sum(young_forest$GX20 > 0 & young_forest$GX20 < 980 & young_forest$GY20 > 0 & young_forest$GY20 < 480) * 400

notyoung_area <- 5e5 - young_area # 48.08 ha
notyoung_area_core <- 5e5 - edge_area - young_area_core # 42.84 ha
core_area <- 5e5 - edge_area # 44.16 ha

# Save all the edited files (this intermediate file is not saved in git repo)
save(list = c('growth8590','growth9095',grep('bci.', ls(), value = TRUE)), file = 'data/bciqcrun.RData')
