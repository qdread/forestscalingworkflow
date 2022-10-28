# Script to draw all plots in manuscript, including main text and supplemental figures

# Specify which model fits to use (3-segment density and 1-segment production)
DENS <- 3
PROD <- 1
year_to_plot <- 1995 # Use results from 1995.

# Load packages and set plotting options ----------------------------------

library(forestscaling) 
library(tidyverse)
library(scales)
library(RColorBrewer)
library(grid)
library(gtable)
library(reshape2)
library(hexbin)
library(Hmisc, pos = 100)

guild_fills <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_fills2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

fg_labels <- c("fg1", "fg2","fg3", "fg4", "fg5")
fg_labels2 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5")
fg_labels3 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5", "unclassified")

guild_names <- paste('fg', 1:5, sep = '')
guild_labels <- c('Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels2 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels3 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium',  "Unclassified")


guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                        fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))

geom_size <- 4


guide <- guides(fill = guide_legend(title = NULL), override.aes = list(fill = NA), color = F)
guide2 <- guides(color = guide_legend(title = NULL))

#---------------------------- Data ---------------------------- 
load('data/rawdataobj1995.RData')
fp <- 'data/data_forplotting'

for (i in dir(fp, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}
################################################################################################
# ------------------------------ Fig 1: Light Interception ---------------------------------
################################################################################################

lightperareacloudbin_fg <- read.csv(file.path(fp, 'lightperareacloudbin_fg.csv'), stringsAsFactors = FALSE)
lightpervolcloudbin_fg <- read.csv(file.path(fp, 'lightpervolcloudbin_fg.csv'), stringsAsFactors = FALSE)
unscaledlightbydbhcloudbin_fg <- read.csv(file.path(fp, 'unscaledlightbydbhcloudbin_fg.csv'), stringsAsFactors = FALSE)

exl <- expression(atop('Light per Crown Area', paste('(W m'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W m'^-3, ')')))
exd <- 'Diameter (cm)'


#----------------------   Fig 1a: Light per crown volume by diameter -----------------------------

ggplot() + 
  geom_point(alpha = 0.01, data = alltree_light_95, 
             aes(x = dbh_corr, y = light_received_byarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareacloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_plant_small(legend = TRUE) 


#----------------------   Fig 1b: Light per crown volume by diameter -----------------------------

 ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, 
             aes(x = dbh_corr, y = light_received_byvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(1, 10, 100)) +
  theme_plant_small(legend = TRUE) 

################################################################################################
# ------------------------------ Fig 3: Life Histories ---------------------------------
################################################################################################

#### Fig 3a, PCA

fgbci <- read_csv('data/functionalgroups.csv')


Fig_3a <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = geom_size, color = "black") + 
  labs(x = 'Slow - Fast Tradeoff', y = 'Recruitment - Stature Tradeoff') +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors)+
  scale_fill_manual(values = guild_fills, labels = guild_labels) + 
  theme_plant_small(legend = TRUE) + guide
Fig_3a

###Fig 3b Growth by Light

title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(atop('Growth per Crown Area', paste('(kg yr'^-1, ' m'^-2,')', sep='')))

obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv('data/data_piecewisefits/lightbyarea_paramci_by_fg.csv', stringsAsFactors = FALSE)

obs_limits <- obs_light_binned %>%
  group_by(fg, year) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

pred_light <- pred_light %>%
  left_join(obs_limits) %>%
  filter(light_area >= min_obs & light_area <= max_obs)

pred_light_5groups <- pred_light %>% filter(!fg %in% c('alltree','unclassified'))

melt_pars <- melt(param_ci, id.vars=1:3)
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

dodge_width <- 0.03
error_bar_width <- 0.04

# correct the error bar width for the number of groups in each bin
obs_light_binned_plotdata <- obs_light_binned %>% 
  filter(year == year_to_plot, mean_n_individuals >= 20, !fg %in% c('alltree', 'unclassified')) %>%
  arrange(factor(fg, levels = c('all', 'fg5','fg4','fg2', 'fg3', 'fg1'))) %>%
  group_by(bin_midpoint, year) %>% 
  mutate(width = sum(c('fg1','fg2','fg3','fg4','fg5') %in% fg)) %>% 
  ungroup

Fig_3b <- ggplot(obs_light_binned_plotdata) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, fill = fg), 
              alpha = 0.4, show.legend = F) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot),
            aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, 
                    group = fg, color = fg, width = 0),
                position = position_dodge(width = dodge_width)) + 
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg),
             size = 4, shape = 21, position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, position = "right", breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c( 0.01, 0.03, 0.1, 0.3)) +
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_fills, labels = guild_labels) +
  theme_plant_small(legend = TRUE)  +guide
Fig_3b

###Fig 3c, Mortality by Light 

bin_mort <- read_csv(file.path(fp, 'obs_mortalitybins.csv')) # Load binned data
bin_mort <- bin_mort %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg2', 'fg3', 'fg1'))) 
mort <- read_csv(file.path(fp, 'obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv('data/data_piecewisefits/mortality_ci_by_fg.csv') # Load fitted values

# Truncate fitted mortality lines to not exceed the range of the observed data, using 20 individuals as the cutoff.
obs_range_mort <- bin_mort %>% 
  filter(variable %in% 'light_per_area', lived + died >= 20) %>%
  group_by(fg) %>%
  summarise(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

fitted_mort_trunc <- fitted_mort %>%  arrange(factor(fg, levels = c('all', 'fg5','fg4','fg2', 'fg3', 'fg1'))) %>%
  left_join(obs_range_mort) %>%
  filter(light_per_area >= min_obs & light_per_area <= max_obs)


Fig_3c <- ggplot(data = fitted_mort_trunc %>% 
                   mutate(fg = factor(fg, labels = guild_labels))) +
  geom_ribbon(aes(x = light_per_area, ymin = q025, ymax = q975, group = fg, fill = fg), 
              alpha = 0.3, show.legend = F) +
  geom_line(aes(x = light_per_area, y = q50, group = fg, color = fg)) +
  geom_point(data = bin_mort %>% 
               filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), (lived+died) >= 20)  %>% 
               mutate(fg = factor(fg, labels = guild_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg),
             shape = 21, size = geom_size) +
  scale_x_log10(name = parse(text = 'Light~per~Crown~Area~(W~m^-2)'), breaks = c(3, 30, 300), limits = c(1.5, 412)) +
  scale_y_continuous(trans = "logit", position = "right", breaks = c(0.03, 0.1, 0.3, 0.6), 
                     labels = c(0.03, 0.1, 0.3, 0.6), limits = c(0.02, 0.65),
                     name = expression(paste("Mortality (5 yr"^-1,")"))) +
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_fills) +
  theme_plant_small(legend = TRUE)  + guide
Fig_3c


########################################################################################
# ------------------------------- Fig 4 Scaling Plots ------------------------------------
########################################################################################

#---- Fig 4a: Growth scaling

Fig_4a <- plot_prod(year_to_plot = 1995,
                fg_names =  c("fg5", "fg4", "fg2", "fg3", "fg1"),
                model_fit = PROD,
                x_limits = c(1, 230),
                y_limits = c(0.001, 2000),
                y_breaks = c(0.001,0.1, 10, 1000),
                plot_errorbar = T,
                error_min = 'q25',
                error_max = 'q75',
                error_bar_width = 0,
                y_labels = c(0.001,0.1,10,1000),
                dodge_width = 0.05)
Fig_4a + theme_plant_small(legend = TRUE)  + guide + 
  theme(axis.text.x = element_text(color = "black", size = 15),  
        axis.ticks.x = element_line(color = "black")) + 
  scale_fill_manual(values = guild_fills, labels = guild_labels) 

#---- Fig 4b: Population Density Scaling

Fig_4b <- plot_dens(year_to_plot = 1995,
                fg_names = fg_labels2,
                model_fit = DENS,
                preddat = pred_dens, 
                x_limits = c(.8, 230),
                y_limits = c(0.001, 20000),
                x_breaks = c(1, 10, 100),
                y_labels = c(0.001, 0.1, 10,1000),
                y_breaks = c(0.001, 0.1,  10, 1000))
Fig_4b + theme_plant_small(legend = TRUE)  + guide +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 


#---- Fig 4c: Total growth scaling
obs_totalprod <- obs_totalprod %>%
  filter(bin_count >= 20)
Fig_4c <- plot_totalprod(year_to_plot = 1995,
                     fg_names = fg_labels2,
                     model_fit_density = DENS, 
                     model_fit_production = PROD,
                     x_limits = c(0.7,230),
                     y_limits = c(0.5, 200),
                     y_breaks = c(0.1, 1, 10, 100),
                     y_labels = c(0.1, 1, 10, 100),
                     preddat = fitted_totalprod)
Fig_4c  + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,10,100), limits = c(0.8, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", breaks = c(3, 10, 30))) +
  theme_plant_small(legend = TRUE)  + guide +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 


########################################################################################
# ------------------------------- Fig 5 Symmetry Plots ------------------------------------
########################################################################################


fitted_indivlight <- read.csv(file.path(fp, 'fitted_indivlight.csv'), stringsAsFactors = FALSE)
fitted_totallight <- read.csv(file.path(fp, 'fitted_totallight.csv'), stringsAsFactors = FALSE)
indivlightbins_fg <- read.csv(file.path(fp, 'obs_indivlight.csv'), stringsAsFactors = FALSE)
totallightbins_fg <- read.csv(file.path(fp, 'obs_totallight.csv'), stringsAsFactors = FALSE)

# Fig 5a      
totallightbins_fg <- totallightbins_fg %>%
  filter(bin_count >= 20)
Fig_5a <- plot_totalprod(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                            model_fit_density = DENS, 
                            model_fit_production = PROD,
                            x_limits = c(0.9,150),
                            y_limits = c(100, 200000),
                            geom_size = 4,
                            y_breaks = c(100, 1000, 10000, 100000),
                            y_labels = c("0.1", "1", "10", "100"), #convert from watts to kilowatts
                            y_name = expression(atop('Total Light Intercepted' , paste('(kW cm'^-1,' ha'^-1,')'))),
                            preddat = fitted_totallight,
                            obsdat = totallightbins_fg,
                            plot_abline = FALSE)
Fig_5a + theme_plant_small(legend = TRUE)  + guide + 
  scale_y_continuous(trans = "log10", position = "left",
                  limits= c(100, 200000), breaks = c(100, 1000, 10000, 100000),
                  labels =  c("0.1", "1", "10", "100"),
                  name = expression(atop('Total Light Intercepted' , paste('(kW cm'^-1,' ha'^-1,')')))) +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 


# Fig 5b
params <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == DENS) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)

growth_slopes <- read.csv('data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
light_slopes <- read.csv('data/data_piecewisefits/light_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_incoming_light', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fg_labels3, labels = guild_labels3))

grob_b <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Light Capture", x = 0.65, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.65, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.35, y = 0.51, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) 

# Plot
Fig_5b <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), aes(x = fg, y = q50, ymin = q025, ymax = q975, fill =  variable, color =variable)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_point(position = position_dodge(width = 0.6), shape = 21, size = 4, color = "black", stroke = 0.5) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  labs( x = NULL, y = 'Scaling Slope') +
  scale_y_continuous(limits = c(-1.05, 1.3)) +
  scale_fill_manual(values = c('gold1', 'darkolivegreen3')) +
  scale_color_manual(values = c('gold3', 'darkgreen')) +
  theme_plant() + theme(axis.text.x = element_text(angle = 25, vjust = 0.7, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3)
Fig_5b



########################################################################################
# ------------------------------- Fig 6, Ratio Scaling ------------------------------------
########################################################################################
# ------------------------------- Fig 6 Relative Abundance & Production ------------------------------------

# Load 1995 bin data
load('data/data_binned/bin_object_singleyear.RData')

prod_ratio_diam <- breeder_stats_bydiam_byyear %>% 
  mutate(ID = 'Short-Tall') %>%
  rename(production_ratio = breeder_production_ratio, density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bydiam_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

prod_ratio_light <- breeder_stats_bylight_byyear %>% 
  mutate(ID = 'Short-Tall') %>%
  rename(production_ratio = breeder_production_ratio, density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bylight_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

# Load fitted values
ratio_fitted_diam <- read_csv('data/data_piecewisefits/ratio_fittedvalues.csv')
ratio_fitted_lightarea <- read_csv('data/data_piecewisefits/ratio_fittedvalues_lightarea.csv')

ratio_fitted_lightarea_prod <- ratio_fitted_lightarea %>%
  filter(variable == 'total production') %>%
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Short-Tall'))

ratio_fitted_diam_density <- ratio_fitted_diam %>%
  filter(variable == 'density', (ratio == 'fast:slow' & dbh < 66) | (ratio == 'pioneer:breeder') & dbh < 16) %>% #limits of binned data
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Short-Tall'))

# ---------------------------------   Fig 6a  Production by light ------------------------------------

prod_ratio <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot() + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_ribbon(aes(x = light_area, ymin = q025, ymax = q975, group = ratio, fill = ratio),
              alpha = 0.3, data = ratio_fitted_lightarea_prod %>% filter(light_area > 7),
              show.legend = F) +
  geom_line(aes(x = light_area, y = q50, group = ratio,  color = ratio), 
            data = ratio_fitted_lightarea_prod %>% filter(light_area > 7)) +
  geom_point(aes(x = bin_midpoint, y = production_ratio, fill = ID), 
             shape = 21, size = 4,  stroke = .5, color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
  theme_plant_small(legend = TRUE)  + 
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(2,330), breaks=c(3,  30,  300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = "Production Ratio") + guide
prod_ratio

#----------------------------- Fig 6B Abundance by Diameter ----------------------------
dens_ratio <- prod_ratio_diam %>% 
  filter(density_ratio > 0) %>%
  filter(n_individuals >= 20) %>%
  ggplot() +
  geom_ribbon(aes(x = dbh, ymin = q025, ymax = q975, group = ratio, fill = ratio), 
              alpha = 0.4, data = ratio_fitted_diam_density, show.legend = F) +
  geom_line(aes(x = dbh, y = q50, group = ratio, color = ratio), data = ratio_fitted_diam_density) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point(aes(x = bin_midpoint, y = density_ratio, fill = ID),
             shape = 21, size = 4,  stroke = .5,  color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
  scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
  scale_x_log10(limits = c(1,100), breaks = c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels = signif, breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = expression("Density Ratio")) + 
  theme_plant_small(legend = TRUE)  + guide
dens_ratio 


########################################################################################
# ------------------------------- Supplementals ------------------------------------
########################################################################################

#---------------- Ext. Fig 1: LAI and Transmittance  ----------------
lai_0 <- read_csv('data/LAI_Depth_Kitajima2005.csv')
pfd_0 <- read_csv('data/PFD_LAI_Kitajima2005.csv')
lai <- lai_0 %>% filter(Species != "Cecropia") 
pfd <- pfd_0 %>% filter(Species != "Cecropia")

fill_sp <- c("darkolivegreen4", "cadetblue2",  "brown3", "gray22") 
color_sp <- c("darkolivegreen", "cadetblue3",  "brown4", "black") 

# Ext. Fig 1a: LAI vs Crown Depth
depth <- ggplot(data = lai, 
                aes(x = Depth, y = LAI, fill = Species, fill = Species, color = Species )) +
  geom_smooth(method = "lm", alpha = 0.3, 
              aes(fill = Species, color = Species),
              show.legend = FALSE) +
  geom_point(aes(fill = Species), size = 4, shape = 21, color = "black") + 
  theme_plant_small(legend = TRUE)  +
  scale_y_log10(limits = c(0.6, 10),labels = signif, breaks = c(1, 3, 10)) +
  scale_x_log10(limits = c(0.2, 15), breaks = c( 0.3, 1, 3, 10), 
                labels = signif, name = "Crown Depth (m)") +
  scale_color_manual(values = color_sp) +
  scale_fill_manual(values = fill_sp) + guide 
depth

# Ext. Fig 1b: Transmittance vs LAI
trans <- ggplot(data = pfd,
                aes( x = LAI, y = PFD, fill = Species, color = Species )) +
  geom_smooth(method = "lm", alpha = 0.3, 
              aes(fill = Species, color = Species), show.legend = F) +
  geom_point(aes(fill = Species), size = 4, shape = 21, color = "black") + 
  theme_plant_small(legend = TRUE)  +
  scale_y_log10(limits = c(2.5, 200), name = "% PFD Transmittance") +
  scale_x_continuous(limits = c(-0.3, 8)) +
  scale_color_manual(values = color_sp) +
  scale_fill_manual(values = fill_sp) + guide 

trans
# -------Ext. Fig 3: see code for Fig 2a and Ruger et al 2018, Fig 1a. note x axis flipped so increases from slow to fast  -----------




# ------------------------ Ext. Fig 2: Heat Map of Growth per Crown Area ~ Light per Crown Area -------------------------
#- note: over 13,000 data points means patterns are obscured if raw data is plotted
# Axis titles
title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(atop('Growth per Crown Area', paste('(kg yr'^-1, ' m'^-2,')', sep='')))

fg_labeler <- c("fg1" = "Fast", "fg2" = "Tall", "fg3" = "Slow", "fg4" = "Short", "fg5" = "Medium")
                     
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))


growth_light_hex <- ggplot(obs_light_raw %>% 
                         filter(year == year_to_plot, fg %nin% c('alltree','unclassified'))) +
  facet_wrap(~ fg, ncol = 2, labeller = as_labeller(fg_labeler)) + 
  geom_hex(aes(x = light_area, y = production_area)) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50), size = 1, color = 'black') +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, labels=signif) +
  hex_scale_log_colors + 
  theme_plant_small() +
  theme(strip.text = element_text(size=14))+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=15))
growth_light_hex


# Overlaying binned means

obs_light_binned2 <- obs_light_binned %>% 
  filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))
growth_light_hex2 <- ggplot(obs_light_raw %>% 
                             filter(year == year_to_plot, fg %nin% c('alltree','unclassified'))) +
  facet_wrap(~ fg, ncol = 2, labeller = as_labeller(fg_labeler)) + 
  geom_hex(aes(x = light_area, y = production_area)) +
  geom_line(data = pred_light_5groups %>% 
              filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
            aes(x = light_area, y = q50, group = fg), color = "black" ,size = 0.5) +
  geom_segment(data = obs_light_binned %>% 
                 filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')),
               aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  geom_point(data = obs_light_binned %>% 
               filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')),
                      shape = 21, size = 2, aes(x = bin_midpoint, y = mean)) +
  scale_x_log10(name = title_x,  limits = c(1, 1000), breaks = c(1,10,100, 1000)) + 
  scale_y_log10(name = title_y, limits = c(0.001, 3), labels = signif, 
                breaks = c(0.001, 0.01,0.1,1,10)) +
  hex_scale_log_colors + 
  theme_plant_small() +
  theme(strip.text = element_text(size=14))+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 15))
growth_light_hex2




#-------------------------- Ext. Fig 4a: Growth Responsiveness to light --------------------------

param_ci$fg <- factor(param_ci$fg ,levels = c("fg1", "fg2", "fg5", "fg3", "fg4"))
guild_labels2_ <- c("Fast", "Tall", "Medium", "Slow", "Short")

ggplot(param_ci %>% filter(fg != 'NA', year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
            aes(x = fg, y = q50, ymin = q025, ymax = q975)) + 
  geom_errorbar(width = 0.4) + geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 25,  vjust = 0.7))+
  scale_x_discrete(name = 'Life History Guild', labels = guild_labels2_) +
  scale_y_continuous(expression(atop('Max. Growth Responsiveness',paste('to Light (kg yr'^-1, ' m'^-2,')'))), 
                     limits = c(0.6, 1.1),
                     breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
  theme_plant() + theme(aspect.ratio = 0.75)


#-------------------------- Ext. Fig 4b: Mortality Responsiveness to light ----------------------
mortal <- read_csv('data/data_forplotting/clean_parameters_mortality.csv')
mortal <- mortal %>% filter(fg != "--")

mortal$fg <- factor(mortal$fg ,levels = c("fast", "large pioneer", "medium", "slow", "small breeder"))
guild_labels2_ <- c("Fast", "Tall", "Medium", "Slow", "Short")
mort_slope <- ggplot(mortal %>% filter(parameter %in% 'slope'),
                     aes(x = fg, y = mean, ymin = q025, ymax = q975)) + 
  geom_errorbar(width = 0.4) + geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 25,  vjust = 0.7))+
  scale_x_discrete(name = 'Life History Guild', labels = guild_labels2_) +
  scale_y_continuous(limits = c(-1.3, -0.2), breaks = c(-1.25, -0.75, -0.25), 
                     expression(atop('Mortality Responsiveness to Light', paste('(yr'^-1,' W'^-1,' m'^-2,')')))) +
  theme_plant() + theme(aspect.ratio = 0.75) 
mort_slope


#-------------------------- Ext. Fig 5: Heat Map of Growth ~ Diameter Scaling -----------------------


#--------------------- Growth by LH Hex ----------------------------
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), 
                                function(x, y) cbind(year = y, x %>% filter(!recruit) %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))


hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))

p0 <- plot_prod_withrawdata(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                            full_names = c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'Unclassified'),
                            x_limits = c(0.7, 316),
                            x_breaks = c(1, 10, 100),
                            y_limits = c(0.0007, 1000),
                            y_breaks = c(.001, .1, 10, 1000),
                            line_types = c('solid', 'dashed'),
                            hex_scale = hex_scale_log_colors,
                            plot_abline = FALSE,
                            plot_fits = PROD)

p <- p0 +  theme_plant_small() + theme(legend.position = 'right', legend.text=element_text(size = 13), 
                legend.title = element_text(size = 15)) +
  scale_y_log10(labels = c(0.001, 0.1, 10, 1000), limits = c(0.0007, 1000),
                breaks = c(0.001, 0.1, 10, 1000), name = expression(paste('Growth (kg yr'^-1,')'))) +
  guides(linetype = FALSE)
p


#----------------------- Ext. Fig 6a: total production from 1985 - 2010 --------------------------
minmax_prod_bycensus <- obs_totalprod %>%
  filter(bin_value > 0, bin_count >= 20, !fg %in% 'unclassified') %>%
  group_by(fg, bin_midpoint) %>%
  summarize(range_min = min(bin_value), range_max = max(bin_value))

p <- ggplot(minmax_prod_bycensus, 
            aes(x = bin_midpoint, ymin = range_min, ymax = range_max, color = fg)) +
  geom_errorbar(size = 1) +
  scale_x_log10(name = 'Diameter (cm)', breaks = c(1,3,10,30,100,300)) + 
  scale_y_log10(expression(paste('Production (kg ha'^-1,' cm'^-1,' yr'^-1,')')),
                breaks = 10^(-2:3), labels = as.character(10^(-2:3)), limits = c(0.1, 200)) +
  scale_color_manual(values = guild_colors2, labels = guild_labels2) +
  theme_plant_small(legend = T) + guide2
p

#-------------- Ext. Fig 6b: Total Crown volume Scaling -----------------------------

totalvolbins_fg <- read.csv(file.path(fp, 'obs_totalvol.csv'), stringsAsFactors = FALSE)

totalvolbins_fg <- totalvolbins_fg %>%
  filter(bin_count >= 20)
p <- plot_totalprod(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = DENS, 
                    model_fit_production = PROD,
                    x_limits = c(0.9, 150),
                    y_limits = c(8, 5500),
                    y_breaks = c(1, 10, 100, 1000),
                    y_labels = c(1, 10, 100, 1000),
                    y_name = expression(atop('Total Crown Volume',paste('(m'^3, ' cm'^-1,' ha'^-1,')'))), 
                    preddat = fitted_totalvol,
                    obsdat = totalvolbins_fg, 
                    plot_abline = FALSE,
                    geom_size = 4)
p +theme_plant_small(legend = TRUE)  + guide +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 


#----------------------- Ext. Fig 7 : Individual light capture scaling --------------------------

alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

p <- ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  theme_plant() + theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))# +
p



# --------------------Ext. Fig 8: Mean PCA scores   -----------------------------
##-- Ext. Fig 8a: Mean score  by Light

fastslowscore_bin_bylight <- fastslowscore_bin_bylight_byyear %>%
  left_join(fastslow_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Survivorship-Growth")

breederscore_bin_bylight <- breederscore_bin_bylight_byyear  %>%
  left_join(breeder_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Recruitment-Stature")

PCA_score_by_light <- rbind(fastslowscore_bin_bylight,breederscore_bin_bylight)

error_bar_width <- .15

PCA_light <- PCA_score_by_light %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_errorbar(width = error_bar_width) +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black") +
  scale_fill_manual(values = c("Recruitment-Stature" = "black", "Survivorship-Growth" = "grey")) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") + 
  scale_x_log10(limits=c(1.8,450), breaks = c(3, 30, 300),
                name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = 'PCA Score')  +
  theme_plant_small(legend = T) + guide
PCA_light


#---- Ext. Fig 8b: PCA by Diameter

fastslowscore_bin_bydiam <- fastslowscore_bin_bydiam_byyear %>%
  left_join(fastslow_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Survivorship-Growth")

breederscore_bin_bydiam <- breederscore_bin_bydiam_byyear %>%
  left_join(breeder_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = 'Recruitment-Stature')

score_bin_bydiam <- rbind(fastslowscore_bin_bydiam, breederscore_bin_bydiam)


PCA_diam <- score_bin_bydiam %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width)  +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Recruitment-Stature" = "black", "Survivorship-Growth" = "grey"))+
  scale_x_log10(name = 'Diameter (cm)', limits = c(.8,100)) + 
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = "PCA Score") +
  theme_plant_small(legend = T) + guide
PCA_diam


#------------------ Ext. Fig 9:Light Scaling---------------------

# Load plotting data
load('data/data_forplotting/light_scaling_plotting_data.RData')
fg_names <- c("Fast", "Tall", "Slow", "Short")
fill_scale <- scale_fill_manual(values = guild_fills[1:4], name = NULL, labels = fg_names, guide = guide_legend(override.aes = list(shape = 21)))
color_scale <- scale_color_manual(values = guild_colors[1:4], name = NULL, labels = fg_names, guide = FALSE)

area_core <- 42.84

obs_indivprod_lightarea$fg <- factor(obs_indivprod_lightarea$fg , labels = c("Fast", "Tall", "Slow", "Short"))
#------------ Ext. Fig 9a: Growth with light
l_growth <- ggplot() +
  geom_ribbon(data = prod_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), 
              alpha = 0.3, show.legend = F) +
  geom_line(data = prod_pred_dat_lightarea %>% filter(light_area >= 7), 
            aes(x = light_area, y = q50, group = fg, color = fg), show.legend = F) +
  geom_point(data = obs_indivprod_lightarea %>% filter(mean_n_individuals >= 20), 
             aes(x = bin_midpoint, y = median, group = fg, fill = fg), 
             shape = 21, color = 'black', size = 4, show.legend = T) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)')) +
  scale_y_log10(labels = signif, limits = c(0.01, 100), name = parse(text = 'Growth~(kg~y^-1)')) +
  fill_scale +
  color_scale +
  theme_plant_small(legend = T)
l_growth

#------------ Ext. Fig 9b: Abundance with light
l_abun <- ggplot() +
  geom_ribbon(data = dens_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), 
              alpha = 0.3, show.legend =F ) +
  geom_line(data = dens_pred_dat_lightarea %>% filter(light_area >= 7), 
            aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_dens_lightarea %>% filter(bin_count >= 20, bin_value > 0), 
             aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg),
             shape = 21, color = 'black', size = 4) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)')) +
  scale_y_log10(labels = signif, limits = c(0.01, 200), name = parse(text = 'Abundance~(ha^-1~cm^-1)')) +
  theme_plant_small(legend = T) +
  fill_scale +
  color_scale
l_abun 

#------------ Ext. Fig 9c: Production with light
l_prod <- ggplot() +
  geom_ribbon(data = totalprod_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), show.legend = F, alpha = 0.3) +
  geom_line(data = totalprod_pred_dat_lightarea %>% filter(light_area >= 7), 
            aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_totalprod_lightarea %>% filter(bin_count >= 20, bin_value > 0), 
             aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg),
             shape = 21, color = 'black', size = 4, show.legend = T) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)')) +
  scale_y_log10(labels = signif, limits = c(0.01, 10), 
                name  = expression(atop('Production', paste('(kg yr'^-1,' cm'^-1,' ha'^-1,')')))) +
  theme_plant_small(legend = T) +
  fill_scale + 
  color_scale

l_prod 

#------------ Ext. Fig 9d: Production ratio with light - See Fig 6a


