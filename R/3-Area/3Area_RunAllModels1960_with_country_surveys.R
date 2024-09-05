#'
#' Run all 3A models
#'

# source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)

## Age-frequencey is direct ageing estimator vs Age-length key estimator
dat_path = here("Data", "3-Area")
fig_path = here("Output", "3-Area")
if(!dir.exists(fig_path)) dir.create(fig_path)

## Age-frequencey is direct ageing estimator vs Age-length key estimator
AF_direct_ageing = T
redo_sim_obs = F
max_N_eff = 500 ## during the initial input calculations some of the input sample sizes were huge, which I don't want effecting the model
min_N_eff = 50
N_eff_multiplier = 0.4

## read in observational datasets
fixed_gear_AF_alk_pooled = readRDS(file = here(dat_path, "Observer_ALK_AF_w_eff.RDS"))
fixed_gear_AF_direct = readRDS(file = here(dat_path, "Observer_direct_AF_w_eff.RDS"))
trawl_observer_LF = readRDS(file = here(dat_path, "Observer_trawl_LF_w_eff.RDS"))
fixed_observer_LF = readRDS(file = here(dat_path, "Observer_fixed_LF_w_eff.RDS"))
survery_AF_direct = readRDS(file = here(dat_path, "Survey_direct_AF_w_eff.RDS"))  %>% mutate(country = ifelse(year <= 1993, "Japan", "United States"))
survery_AF_alk_pooled = readRDS(file = here(dat_path, "Survey_ALK_AF_w_eff.RDS"))  %>% mutate(country = ifelse(year <= 1993, "Japan", "United States"))

fixed_gear_AF_direct = fixed_gear_AF_direct %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
fixed_gear_AF_alk_pooled = fixed_gear_AF_alk_pooled %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
survery_AF_direct = survery_AF_direct %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
survery_AF_alk_pooled = survery_AF_alk_pooled %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
trawl_observer_LF = trawl_observer_LF %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
fixed_observer_LF = fixed_observer_LF %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))

fixed_gear_AF_direct = fixed_gear_AF_direct %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
fixed_gear_AF_alk_pooled = fixed_gear_AF_alk_pooled %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
survery_AF_direct = survery_AF_direct %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
survery_AF_alk_pooled = survery_AF_alk_pooled %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
trawl_observer_LF = trawl_observer_LF %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
fixed_observer_LF = fixed_observer_LF %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
design_survey_index = readRDS(file = file.path("Data", "Survey", "regional_abundance_estimates.RDS"))
design_survey_index$area_lab = design_survey_index$NPFMC.Sablefish.Management.Area

design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "Aleutians" ~ "BS_AI_WGOA",
                                                                 area_lab == "Bering Sea" ~ "BS_AI_WGOA",
                                                                 area_lab == "Western Gulf of Alaska" ~ "BS_AI_WGOA",
                                                                 area_lab == "Central Gulf of Alaska" ~ "CGOA",
                                                                 area_lab == "Eastern Gulf of Alaska" ~ "EGOA",
                                                                 TRUE ~ area_lab))
## sum estimates for BS, AI, WGOA
design_survey_index = design_survey_index %>% group_by(Country, Year, area_lab) %>% summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)
design_survey_index$type = "Updated"
ggplot(design_survey_index, aes(Year, sum_estimates, col = Country, linetype = Country, fill = Country, alpha = 0.3)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI)) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Relative index (RPN)", linetype = "", col = "", fill = "") +
  guides(alpha = "none") +
  facet_wrap(~area_lab, scales = "free_y") +
  theme_bw()


## Tag data
tag_recovery_df = readRDS(file = file.path(dat_path, "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = file.path(dat_path, "Tag_release_summarised.RDS"))

## read in Catch
full_catch_df = readRDS(file = file.path(dat_path, "Catch_by_year_area_gear.RDS"))

## read in historic extrapolated Catch
fixed_gear_with_imputation = readRDS(file = here("Data", "3-Area", "fixed_gear_with_imputations_S1.RDS"))
trawl_gear_with_imputation = readRDS(file = here("Data", "3-Area", "trawl_gear_with_imputations_S1.RDS"))

## bring ADMB inputs to help configure initial model
sab_curr <- dget(here("Data", "ADMB", "tem.rdat")) 
sab_rep <-readLines(here("Data", "ADMB", "sable.rep"))
sab_inputs <- digest_sab_input(sablefish_input_filename = here("Data", "ADMB", "tem_2022_na_wh.dat")) 
sab_par <- readLines(here("Data", "ADMB", "tem.par")) 
sab_ctl <- digest_sab_ctl(sablefish_control_filename = here("Data", "ADMB", "tem.ctl")) 
## Model dimensions

years = 1960:2021
yrs_to_add = 1977 - 1960

min_age = 2
max_age = 31
ages = min_age:max_age
region_key = data.frame(area = c("BS_AI_WGOA","CGOA","EGOA"), TMB_ndx = c(0:2))
length_bins = seq(from = 41, to = 99,by = 2)
M = 0.104884
reg_lvls = region_key$area[region_key$TMB_ndx + 1] ## need to be in correct order
year_lvls = years ## need to be in correct order
## length bins
sex_length_lvls = paste0(rep(c("M","F"), each = length(length_bins)), length_bins)
sex_age_lvls = paste0(rep(c("M","F"), each = length(ages)), ages)
trawl_observer_LF$sex_length = factor(trawl_observer_LF$sex_length, levels = sex_length_lvls, ordered = T)
fixed_observer_LF$sex_length = factor(fixed_observer_LF$sex_length, levels = sex_length_lvls, ordered = T)
fixed_gear_AF_alk_pooled$sex_age = factor(fixed_gear_AF_alk_pooled$sex_age, levels = sex_age_lvls, ordered = T)
fixed_gear_AF_direct$sex_age = factor(fixed_gear_AF_direct$sex_age, levels = sex_age_lvls, ordered = T)
survery_AF_direct$sex_age = factor(survery_AF_direct$sex_age, levels = sex_age_lvls, ordered = T)
survery_AF_alk_pooled$sex_age = factor(survery_AF_alk_pooled$sex_age, levels = sex_age_lvls, ordered = T)
#'
#' TMB data
#'
data <- list()
data$model = "TagIntegrated"
data$ages = ages
data$years = years
data$length_bins = length_bins
data$n_regions = nrow(region_key)
data$n_movement_time_blocks = 1
n_regions = data$n_regions
survey_labels = c("Japanese", "Domestic")
data$n_surveys = length(survey_labels)
n_ages = length(data$ages) 
n_length_bins = length(data$length_bins) # the last length bin value is the minimum for a length plus group
data$n_projections_years = 0
data$do_projection = 0
n_projyears = length(data$years) +  data$n_projections_years
n_years = length(data$years)
projyears = min(data$years):(max(data$years) + data$n_projections_years)

data$global_rec_devs = 1
data$rec_devs_sum_to_zero = 0
data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(n_years);
data$n_init_rec_devs = 0 ## number of non-equilibrium age-structure
data$M = matrix(0.104884, nrow = n_ages, ncol = n_projyears)
data$maturity = matrix(sab_curr$growthmat[,5], nrow = n_ages, ncol = n_projyears, byrow = F)

## weight at age
weight_at_age_f <-as.numeric(unlist(strsplit(sab_rep[(grep("Weight Females",sab_rep)+1):(grep("Weight Females",sab_rep)+n_years)],split=" "))) #Get numbers at age
weight_at_age_m<-as.numeric(unlist(strsplit(sab_rep[(grep("Weight Males",sab_rep)+1):(grep("Weight Males",sab_rep)+n_years)],split=" "))) #Get numbers at age
weight_at_age_f = subset(weight_at_age_f, subset = !is.na(weight_at_age_f))
weight_at_age_m = subset(weight_at_age_m, subset = !is.na(weight_at_age_m))

# turn to matrix
weight_at_age_f_mat = matrix(weight_at_age_f, byrow = T, nrow = n_years, ncol = n_ages)
weight_at_age_m_mat = matrix(weight_at_age_m, byrow = T, nrow = n_years, ncol = n_ages)
##
data$female_mean_weight_by_age = t(weight_at_age_f_mat)
data$male_mean_weight_by_age = t(weight_at_age_m_mat)

## age length transition matrix
data$female_age_length_transition = data$male_age_length_transition = array(0, dim = c(n_ages, n_length_bins, n_years), dimnames = list(data$ages, data$length_bins, data$years))
data$male_age_length_transition[,,as.character(1977:1994)] = sab_inputs$male_length_at_age_60_96
data$female_age_length_transition[,,as.character(1977:1994)] = sab_inputs$female_length_at_age_60_96
data$male_age_length_transition[,,as.character(1995:2020)] = sab_inputs$male_length_at_age_97_22
data$female_age_length_transition[,,as.character(1995:2020)] = sab_inputs$female_length_at_age_97_22


## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
data$fixed_movement_matrix = array(0, dim = c(n_regions, n_regions, 1));
data$movement_matrix = array(0, dim = c(n_regions, n_regions, 1));
data$movement_time_block_indicator = rep(0, n_years)
movement_matrix = fixed_movement_matrix = matrix(0, nrow = n_regions, ncol = n_regions);
diag(fixed_movement_matrix) = 1
diag(movement_matrix) = 0.9
data$fixed_movement_matrix[,,1] = fixed_movement_matrix
movement_matrix = movement_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1)
# renormalise
movement_matrix = sweep(movement_matrix, 1, STATS = rowSums(movement_matrix), "/")
data$movement_matrix[,,1] = movement_matrix

## From Table 1.2 from PhD
fenske_move_matrix = matrix(c(67.5, 22.3 , 10.2,
                              23.9, 37.1, 39.0 ,
                              7.9, 28.2, 63.9), nrow = 3, byrow = T)


data$spawning_time_proportion = rep(0, n_projyears)
data$sigma_R = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("sigr:",sab_par)+1],split=" "))))
data$SrType = 3
data$apply_fixed_movement = 0; ##
data$do_recruits_move = 0

#'
#' Fishing inputs
#' 
data$F_method = 1
data$F_max = 3
data$F_iterations = 4
## 
data$prop_F_hist = 0
# drop years outside model years
full_catch_df = full_catch_df %>% filter(year %in% years)
full_catch_df$year_f = factor(full_catch_df$year, levels = year_lvls, ordered = T)
full_catch_df$region_f = factor(full_catch_df$area, levels = reg_lvls, ordered = T)
## convert catch to metric tonnes
full_catch_df$catch_mt = full_catch_df$catch / 1000
trwl_catch = full_catch_df %>% filter(fmp_gear == "TRW") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)
fixed_catch = full_catch_df %>% filter(fmp_gear == "HAL") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)

## replace 0 values with low values
fixed_catch_real = fixed_catch %>% dplyr::select(-region_f) %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))
trwl_catch_real = trwl_catch %>% dplyr::select(-region_f) %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))

fixed_gear_with_imputation$year_f = factor(fixed_gear_with_imputation$year, levels = data$years, ordered = T)
fixed_gear_with_imputation$region_f = factor(fixed_gear_with_imputation$area, levels = region_key$area[region_key$TMB_ndx + 1], ordered = T)
trawl_gear_with_imputation$year_f = factor(trawl_gear_with_imputation$year, levels = data$years, ordered = T)
trawl_gear_with_imputation$region_f = factor(trawl_gear_with_imputation$area, levels = region_key$area[region_key$TMB_ndx + 1], ordered = T)

trwl_catch = trawl_gear_with_imputation %>% filter(year < 1977) %>% arrange(region_f) %>% ungroup() %>% dplyr::select(imputed_catch, year_f, region_f) %>% pivot_wider(names_from = year_f, values_from = imputed_catch)
fixed_catch = fixed_gear_with_imputation  %>% filter(year < 1977)  %>% arrange(region_f)  %>% ungroup() %>% dplyr::select(imputed_catch, year_f, region_f) %>% pivot_wider(names_from = year_f, values_from = imputed_catch)

trwl_catch_temp = fixed_catch_temp = matrix(0, nrow = data$n_regions, ncol = n_projyears)

trwl_catch_temp[,1:yrs_to_add] = as.matrix(trwl_catch %>% select(!region_f))
trwl_catch_temp[,(yrs_to_add+1):n_projyears] = as.matrix(trwl_catch_real[,(yrs_to_add+1):n_projyears])
fixed_catch_temp[,1:yrs_to_add] = as.matrix(fixed_catch %>% select(!region_f))
fixed_catch_temp[,(yrs_to_add+1):n_projyears] = as.matrix(fixed_catch_real[,(yrs_to_add+1):n_projyears])

trwl_catch_temp[trwl_catch_temp == 0] = 0.001
fixed_catch_temp[fixed_catch_temp == 0] = 0.001

data$fixed_fishery_catch = fixed_catch_temp
data$trwl_fishery_catch = trwl_catch_temp

## control variables
# 2 time-blocks for fixed gear selectivity
data$fixed_sel_type = as.vector(rep(0, 2), mode = "integer")
data$fixed_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")
#data$fixed_sel_by_year_indicator[projyears %in% 1995:2015] = 1
data$fixed_sel_by_year_indicator[projyears > 2015] = 1

# single time-block for trawl fishery
data$trwl_sel_type = as.vector(rep(1, 1), mode = "integer")
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

## two time-blocks for survey selectivity
## three time-blocks for survey selectivity
data$srv_sel_type = matrix(0, nrow = 1, ncol = data$n_surveys)
data$srv_sel_by_year_indicator = matrix(0, nrow = n_projyears, ncol = data$n_surveys)
#data$srv_sel_by_year_indicator[projyears %in% 1988:1994, 1] = 1
#data$srv_sel_by_year_indicator[projyears > 1994, 1] = 2

#'
#' Tag release stuff
#'
tag_release_years = c(1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016)
data$tag_release_event_this_year = rep(0, n_years) ## no tag releases
data$tag_release_event_this_year[data$years %in% tag_release_years] = 1
data$male_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))
data$female_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))
for(y_ndx in 1:length(tag_release_years)) {
  this_release_year = tag_release_df %>% filter(release_year == tag_release_years[y_ndx])
  if(nrow(this_release_year) > 0) {
    regions_to_release = unique(this_release_year$region_release) 
    for(r_ndx in 1:length(regions_to_release)) {
      region_ndx = region_key$TMB_ndx[which(region_key$area %in% regions_to_release[r_ndx])] + 1
      data$male_tagged_cohorts_by_age[, region_ndx, y_ndx] = (this_release_year %>% filter(sex == "M", region_release == regions_to_release[r_ndx]))$Nage_at_release
      data$female_tagged_cohorts_by_age[, region_ndx, y_ndx] = (this_release_year %>% filter(sex == "F", region_release == regions_to_release[r_ndx]))$Nage_at_release
    }
  } else {
    ## no tag release event
    data$tag_release_event_this_year[data$years %in% tag_release_years[y_ndx]] = 0
  }
  
}

#paste(data$years[which(data$tag_release_event_this_year == 1)], collapse = ", ")

data$n_years_to_retain_tagged_cohorts_for = 15
data$initial_tag_induced_mortality = rep(0.1, length(tag_release_years))
data$annual_tag_shedding_rate = 0.02

#'
#' Observation data
#' 
data$ageing_error_matrix = sab_inputs$ageing_error
#data$ageing_error_matrix = matrix(0, nrow = n_ages, ncol = n_ages)
#diag(data$ageing_error_matrix) = 1.0

# reformat year and region into factos so we can keep ordering
if(AF_direct_ageing) {
  fixed_gear_AF = fixed_gear_AF_direct %>% filter(year %in% years)
} else {
  fixed_gear_AF = fixed_gear_AF_alk_pooled %>% filter(year %in% years)
}
fixed_gear_AF$region_f = factor(fixed_gear_AF$region, levels = reg_lvls, ordered = T)
fixed_gear_AF$year_f = factor(fixed_gear_AF$year, levels = year_lvls, ordered = T)

trawl_observer_LF = trawl_observer_LF %>% filter(year %in% years)
trawl_observer_LF$region_f = factor(trawl_observer_LF$region, levels = reg_lvls, ordered = T)
trawl_observer_LF$year_f = factor(trawl_observer_LF$year, levels = year_lvls, ordered = T)
fixed_observer_LF = fixed_observer_LF %>% filter(year %in% years)
fixed_observer_LF$region_f = factor(fixed_observer_LF$region, levels = reg_lvls, ordered = T)
fixed_observer_LF$year_f = factor(fixed_observer_LF$year, levels = year_lvls, ordered = T)

if(AF_direct_ageing) {
  survery_AF = survery_AF_direct %>% filter(year %in% years)
} else {
  survery_AF = survery_AF_alk_pooled %>% filter(year %in% years)
}
survery_AF$region_f = factor(survery_AF$region, levels = reg_lvls, ordered = T)
survery_AF$year_f = factor(survery_AF$year, levels = year_lvls, ordered = T)

## drop years outside of model years. can cause NA's
survey_index = design_survey_index %>% filter(Year %in% years)
survey_index$region_f = factor(survey_index$area_lab, levels = reg_lvls, ordered = T)
survey_index$year_f = factor(survey_index$Year, levels = year_lvls, ordered = T)

indicator_fun = function(x) {
  ifelse(is.na(sum(x)), 0, 1)
}

fixed_AF_indicator = fixed_gear_AF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
trawl_LF_indicator = trawl_observer_LF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
fixed_LF_indicator = fixed_observer_LF %>% filter(!year %in% unique(fixed_gear_AF$year)) %>% ungroup() %>% dplyr::select(region_f, year_f, P) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(P), values_fn = indicator_fun)
survey_AF_indicator = survery_AF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
survey_index_indicator = survey_index %>% ungroup() %>% dplyr::select(region_f, year_f, Year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(Year), values_fn = indicator_fun)

### Fixed gear fishery AF
data$fixed_catchatage_indicator = as.matrix(fixed_AF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$obs_fixed_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years))
## fill the container
obs_years = unique(fixed_gear_AF$year)
obs_reg = unique(fixed_gear_AF$region)
N_eff = 100 ## sample size
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (fixed_gear_AF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx]))
    if(nrow(this_df) > 0)
      data$obs_fixed_catchatage[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N
  }
}

data$fixed_catchatage_covar_structure = 0
data$fixed_catchatage_comp_likelihood = 0

### Trawl gear fishery LF
data$trwl_catchatlgth_indicator = as.matrix(trawl_LF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$obs_trwl_catchatlgth = array(0, dim = c(n_length_bins * 2, n_regions, n_years))
## fill the container
obs_years = sort(unique(trawl_observer_LF$year))
obs_reg = unique(trawl_observer_LF$region)
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = trawl_observer_LF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx])
    if(nrow(this_df) > 0) {
      this_n_eff = unique(this_df$eff_N)
      this_df = this_df %>% dplyr::select(P, sex_length, eff_N) %>% complete(sex_length, fill = list(eff_N =this_n_eff )) %>% arrange(sex_length)
      this_df$P = this_df$P %>% replace_na(0)
      if(sum(this_df$P) == 0) {
        data$trwl_catchatlgth_indicator[this_region_ndx, this_year_ndx] = 0
      } else {
        data$obs_trwl_catchatlgth[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N
      }
    }
  }
}

data$trwl_catchatlgth_covar_structure = 0
data$trwl_catchatlgth_comp_likelihood = 0

### Fixed gear fishery LF
data$fixed_catchatlgth_indicator = array(0, dim = c(n_regions, n_years))
data$obs_fixed_catchatlgth = array(0, dim = c(n_length_bins * 2, n_regions, n_years))

## fill the container
include_fixed_LF = T
obs_years = sort(unique(fixed_observer_LF$year))
obs_reg = unique(fixed_observer_LF$region)
if(include_fixed_LF) { # exclude LF observation
  data$fixed_catchatlgth_indicator = as.matrix(fixed_LF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
  obs_years = unique(fixed_observer_LF$year) 
  
  ## cut years that already in the age comp don't want to double dip
  obs_years = obs_years[!obs_years %in% unique(fixed_gear_AF_alk_pooled$year)]
  obs_reg = unique(fixed_observer_LF$region)
  for(y_ndx in 1:length(obs_years)) {
    this_year_ndx = which(years %in% obs_years[y_ndx])
    for(r_ndx in 1:length(obs_reg)) {
      this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
      this_df = fixed_observer_LF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx])
      if(nrow(this_df) > 0) {
        this_n_eff = unique(this_df$eff_N)
        this_df = this_df %>% dplyr::select(P, sex_length, eff_N) %>% complete(sex_length, fill = list(eff_N =this_n_eff )) %>% arrange(sex_length)
        this_df$P = this_df$P %>% replace_na(0)
        
        if(sum(this_df$P) == 0) {
          data$fixed_catchatlgth_indicator[this_region_ndx, this_year_ndx] = 0
        } else {
          data$obs_fixed_catchatlgth[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N
        }
      }
    }
  }
}
data$fixed_catchatlgth_covar_structure = 0
data$fixed_catchatlgth_comp_likelihood = 0

### Survey LL proportions at age
data$srv_catchatage_indicator = array(0, dim = c(n_regions, n_years,data$n_surveys))
US_survey_years = unique((survery_AF %>% filter(country == "United States"))$year)
jap_survey_years = unique((survery_AF %>% filter(country == "Japan"))$year)

data$obs_srv_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years, data$n_surveys))
obs_years = unique(survery_AF$year)
obs_reg = unique(survery_AF$region)
## populate container
countries = unique(survery_AF$country)
## populate container
for(c_ndx in 1:length(countries)) { ## go backwards becuase its starts with Japan
  for(y_ndx in 1:length(obs_years)) {
    this_year_ndx = which(years %in% obs_years[y_ndx])
    for(r_ndx in 1:length(obs_reg)) {
      this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
      this_df = (survery_AF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx], country == countries[c_ndx])) %>% arrange(sex_age)
      if(nrow(this_df) > 0) {
        data$obs_srv_catchatage[,this_region_ndx, this_year_ndx, c_ndx] = this_df$P * this_df$eff_N
        data$srv_catchatage_indicator[this_region_ndx, this_year_ndx, c_ndx] = 1
      }
    }
  }
}
data$srv_catchatage_covar_structure = 0
data$srv_catchatage_comp_likelihood = rep(0, data$n_surveys)

### Survey LL index
data$obs_srv_bio = data$obs_srv_se = array(0, dim = c(n_regions, n_years, data$n_surveys))
data$srv_bio_indicator = array(0, dim = c(n_regions, n_years,data$n_surveys))

countries = unique(design_survey_index$Country)
obs_years = unique(design_survey_index$Year)
obs_reg = unique(design_survey_index$area_lab)
likelihood_type = "Lognormal"
## populate container
for(c_ndx in 1:length(countries)) {
  for(y_ndx in 1:length(obs_years)) {
    mod_yr_ndx = which(data$years == obs_years[y_ndx])
    for(r_ndx in 1:length(region_key$area)) {
      this_df = design_survey_index %>% filter(Country == countries[c_ndx], Year == obs_years[y_ndx], area_lab == region_key$area[r_ndx])
      if(nrow(this_df) > 0) {
        data$srv_bio_indicator[r_ndx, mod_yr_ndx, c_ndx] = 1
        data$obs_srv_bio[r_ndx, mod_yr_ndx, c_ndx] = this_df$sum_estimates## model numbers are in 1000's is in kilo tonnes
        #data$obs_srv_se[this_region_ndx, this_year_ndx] = this_df$se
        data$obs_srv_se[r_ndx, mod_yr_ndx, c_ndx] = (this_df$se)  # scale by 1000 as well
      }
    }
  }
}

data$srv_bio_likelihood = rep(1, data$n_surveys)
data$srv_obs_is_abundance = rep(1, data$n_surveys)
data$srv_q_by_year_indicator = matrix(0, nrow = n_years, ncol = data$n_surveys)
## three q time-blocks
#data$srv_q_by_year_indicator[projyears %in% 1988:1994, 1] = 1
#data$srv_q_by_year_indicator[projyears > 1994, 1] = 2
data$srv_q_transformation = rep(0, data$n_surveys)
data$q_is_nuisance = rep(0, data$n_surveys)


###
# Tag-recoveries
# This can take a bit of time to get your head around!
###
include_tag_recoveries = T # tags are turned off using evaluate-tag-likelihood switch
include_zero_tag_recovery_events = T
tag_recovery_years = 1978:2016

data$tag_recovery_indicator_by_year = rep(0, n_years) ## no tag releases
data$obs_tag_recovery = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))
data$tag_recovery_indicator = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))

if(include_tag_recoveries) {
  # drop any recovery years before release years
  #tag_recovery_years = tag_recovery_years[which(tag_recovery_years %in% (tag_release_years + 1))] ## the plus one is because we don't allow a recovery unless after a year at release
  #
  
  data$tag_recovery_indicator_by_year[data$years %in% tag_recovery_years] = 1
  tag_recovery_df$release_event_year = tag_recovery_df$recovery_year - tag_recovery_df$release_year
  for(y_ndx in 1:length(tag_recovery_years)) {
    this_recovery_year = tag_recovery_df %>% filter(recovery_year == tag_recovery_years[y_ndx])
    recovery_regions = region_key$area[order(region_key$TMB_ndx)] ## all regions
    for(r_ndx in 1:length(recovery_regions)) {
      region_ndx = region_key$TMB_ndx[which(region_key$area %in% recovery_regions[r_ndx])] + 1
      ## tags must be at liberty for at least one year
      this_recovery_year_region = this_recovery_year %>% filter(release_year <= tag_recovery_years[y_ndx] - 1,release_year %in% tag_release_years, recovery_region == recovery_regions[r_ndx])
      # now find how many 'release-event' tags we recovered during this recovery year and region.- which is release year x release region specific blahh!!!
      for(release_event_ndx in 3:(data$n_years_to_retain_tagged_cohorts_for)) { # index starts at two because 1 indicates they were released this year and we don't consider recoveries until after a year at liberty
        this_release_year = tag_recovery_years[y_ndx] - (release_event_ndx - 1)
        if(!this_release_year %in% tag_release_years)
          next; # there can be early release events in the data that are not in the model which can cause an NaN via 0 expected values in the likelihood calculations
        if(release_event_ndx == (data$n_years_to_retain_tagged_cohorts_for + 1)) {
          this_release_event_df = this_recovery_year_region %>% filter(release_year >= this_release_year) ## we can have a pooled tag group
        } else {
          this_release_event_df = this_recovery_year_region %>% filter(release_year == this_release_year)
        }
        release_regions = region_key$area[order(region_key$TMB_ndx)] ## all regions
        for(rel_ndx in 1:length(release_regions)) {
          ## not all release events in the data are actually released in the model
          ## cannot have a recovery obs because this isn't a release event
          mod_release_ndx = which(tag_release_years %in% this_release_year)
          if((sum(data$male_tagged_cohorts_by_age[,rel_ndx, mod_release_ndx]) + sum(data$female_tagged_cohorts_by_age[,rel_ndx, mod_release_ndx])) <= 0)
            next;
          release_ndx = region_key$TMB_ndx[which(region_key$area %in% release_regions[rel_ndx])] + 1
          model_tag_release_ndx = get_tag_release_ndx(release_ndx, release_event_ndx, n_regions)
          ## there was a possible observation for this release event
          if(include_zero_tag_recovery_events)
            data$tag_recovery_indicator[model_tag_release_ndx, r_ndx, y_ndx] = 1
          ## not all regions were 'release events' check we actually released fish in this event
          specific_recovery_df = this_release_event_df %>% filter(region_release == release_regions[rel_ndx])
          if(nrow(specific_recovery_df) > 0) {
            ## tell the model a recovery observation happened for this release event
            data$tag_recovery_indicator[model_tag_release_ndx, r_ndx, y_ndx] = 1
            ## get the recovery age frequency
            data$obs_tag_recovery[model_tag_release_ndx, r_ndx, y_ndx] = sum((specific_recovery_df %>% group_by(sex_age) %>% summarise(Nage_at_recovery = sum(Nage_at_recovery)))$Nage_at_recovery)
          }
        }
      }
    }
  }
}
## are there observatiosn for the plus group
if(sum(data$tag_recovery_indicator[1:5,,]) != 0)
  cat("found recoveries in the first year of release\n")
if(sum(data$tag_recovery_indicator[(dim(data$tag_recovery_indicator)[1] - data$n_regions + 1):dim(data$tag_recovery_indicator)[1],,]) != 0)
  cat("found recoveries in the plus group\n")

cat("The number of recovery events ", sum(data$tag_recovery_indicator), "\n")

data$tag_likelihood = 1
data$evaluate_tag_likelihood = 0

## Projection period
data$future_recruitment_type = 0
data$year_ndx_for_empirical_resampling = c(0,n_years - 1)
data$future_fishing_type = 0
data$future_fishing_inputs_trwl = array(0.1, dim = c(data$n_regions, data$n_projections_years))
data$future_fishing_inputs_fixed = array(0.1, dim = c(data$n_regions, data$n_projections_years))

#'
#' TMB Parameter definition OM values or starting values for EM
#'
parameters <- list()
parameters$ln_mean_rec = rnorm(data$n_regions, suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_mean_rec:",sab_par)+1],split=" ")))), 0.3)

## Fishery selectivity
parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_trwl_sel_pars = array(0, dim = c(length(unique(data$trwl_sel_by_year_indicator )), 2, 2))
parameters$ln_srv_sel_pars = array(0, dim = c(max(apply(data$srv_sel_by_year_indicator, 2,FUN = function(x){length(unique(x))})), 2, 2, data$n_surveys))

## populate parameters Note some of the male delta values are set to the female values. Line 1800 tem.tpl
if(dim(parameters$ln_fixed_sel_pars)[1] == 3) {
  parameters$ln_fixed_sel_pars[1,1,1]  = parameters$ln_fixed_sel_pars[2,1,1]  = parameters$ln_fixed_sel_pars[3,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = parameters$ln_fixed_sel_pars[3,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = parameters$ln_fixed_sel_pars[3,1,2] =suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] = parameters$ln_fixed_sel_pars[3,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_fixed_sel_pars)[1] == 2) {
  parameters$ln_fixed_sel_pars[1,1,1]  = parameters$ln_fixed_sel_pars[2,1,1]  = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
}
if(F){
  parameters$ln_fixed_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
}
parameters$ln_trwl_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_m:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))

## NOTE: all delta parameters for all survey selectivities are fixed based on srv_1
if(dim(parameters$ln_srv_sel_pars)[1] == 1) {
  parameters$ln_srv_sel_pars[1,1,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1,1:data$n_surveys] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2,1:data$n_surveys] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_sel_pars)[1] == 2) {
  parameters$ln_srv_sel_pars[1,1,1,1:data$n_surveys] =  parameters$ln_srv_sel_pars[2,1,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,2,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,1,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,2,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_sel_pars)[1] == 3) {
  parameters$ln_srv_sel_pars[1,1,1,1:data$n_surveys] =parameters$ln_srv_sel_pars[2,1,1,1:data$n_surveys] = parameters$ln_srv_sel_pars[3,1,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,2,1,1:data$n_surveys] = parameters$ln_srv_sel_pars[3,2,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,1,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[3,1,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[2,2,2,1:data$n_surveys] = parameters$ln_srv_sel_pars[3,2,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
}

## movement pars
parameters$transformed_movement_pars = array(NA, dim = c(n_regions - 1, n_regions, data$n_movement_time_blocks))
for(i in 1:n_regions)
  parameters$transformed_movement_pars[,i,1] = simplex(movement_matrix[i,])

parameters$ln_fixed_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish1:",sab_par)+1],split=" "))))
parameters$ln_fixed_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish1:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_trwl_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish3:",sab_par)+1],split=" "))))
parameters$ln_trwl_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish3:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_init_F_avg = parameters$ln_fixed_F_avg
parameters$trans_srv_q = array(log(7000), dim =c(n_regions, max(apply(data$srv_q_by_year_indicator, 2,FUN = function(x){length(unique(x))})) , data$n_surveys))
admb_ln_rec_devs = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_rec_dev:",sab_par)+1],split=" "))))
admb_ln_rec_devs = subset(admb_ln_rec_devs, subset = !is.na(admb_ln_rec_devs))
data$map_simplex_ycs_estimated = rep(0:(length(data$years) - 1))
data$standardise_ycs = 0
if(data$rec_devs_sum_to_zero == 1) {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
    
  }
} else {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = data$n_regions, ncol =  n_years)
  }
}

#parameters$ln_rec_dev = rep(0 - 0.5 * data$sigma_R * data$sigma_R, n_projyears)
parameters$ln_catch_sd = log(0.02)

if(data$n_init_rec_devs == 0) {
  parameters$ln_init_rec_dev = 0#rev(admb_ln_rec_devs[1:data$n_init_rec_devs]) ## note the reverse fucntion here. see source code for why we have to do this
} else {
  parameters$ln_init_rec_dev = rep(0, data$n_init_rec_devs)
}
parameters$logistic_tag_reporting_rate = matrix(logit(0.50), nrow = data$n_regions, ncol = max(sum(data$tag_recovery_indicator_by_year),1))

parameters$ln_tag_phi = log(0.5)
parameters$ln_sigma_R = log(data$sigma_R)
parameters$ln_sigma_init_devs = log(0.2)

parameters$trans_trwl_catchatlgth_error = log(1)
parameters$trans_fixed_catchatlgth_error = log(1)
parameters$trans_fixed_catchatage_error = log(1)
parameters$trans_srv_catchatage_error = rep(log(1), data$n_surveys)

parameters$logistic_prop_recruit_male = rep(logit(0.5), length(data$years))
parameters$trans_SR_pars = rep(qlogis(0.8), 1)

#'
#' Initial parameter setup that will be shared among model runs
#' these relate to parameters that are passed to the
#' set_up_parameters function, which controls what parameters are estimated
#' or shared. This was done to ensure these were consistent among multiple runs
#' unless a run was specifically investigating/changing this assumption

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_shared_by_sex  = F
fixed_sel_second_shared_by_sex   = T
trwl_sel_first_shared_by_sex  = F
trwl_sel_second_shared_by_sex  = F
recruit_dev_years_not_to_estimate = NULL  ## don't estimate the last one
srv_q_spatial = F
tag_reporting_rate = "constant"
est_init_F = F
est_catch_sd = F
est_movement = T
est_prop_male_recruit = "off"


######################
## Run initial model with
## no tag-recovery data
## this is done by setting evaluate_tag_likelihood = 0
## don't estimate movement. When you try estimate movement you get non PD hessian
## if tagging data is excluded
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_NoTagData_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

# turn off tag data
data$evaluate_tag_likelihood = 0
# are we estimating movement - I have tried to estimate movement with no tagging data but I couldn't 
# get PD covariance matrix
data$apply_fixed_movement = 1


## check parameters
validate_input_data_and_parameters(data, parameters)

## save initial data and parameters
saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))

## plot some input info
input_obs = plot_input_observations(data, region_key = region_key, survey_labels = survey_labels)
# tmp_obs = get_input_observations(data, region_key = region_key)
input_obs
ggsave(filename = file.path(fig_path, "Observation_Frequency.png"), width = 10, height = 10)

time_blocks = plot_input_timeblocks(data)
time_blocks
ggsave(filename = file.path(fig_path, "time_blocks.png"), width = 7, height = 6)

plot_input_catches(data, region_key = region_key)
ggsave(filename = file.path(fig_path, "InputCatches.png"), width = 8, height = 6)


## frequency of tag release and recoveries
BS_AI_WGOA_releases = plot_frequency_of_tag_release_and_recoveries(data, region_key = region_key, release_ndx_to_plot = 1:90, release_region_to_plt = "BS_AI_WGOA")
BS_AI_WGOA_releases
ggsave(filename = file.path(fig_path, "BS_AI_WGOA_tag_release_and_recovery.png"), width = 7, height = 7)
CGOA_releases = plot_frequency_of_tag_release_and_recoveries(data, region_key = region_key, release_ndx_to_plot = 1:90, release_region_to_plt = "CGOA")
CGOA_releases
ggsave(filename = file.path(fig_path, "CGOA_tag_release_and_recovery.png"), width = 12, height = 8)
EGOA_releases = plot_frequency_of_tag_release_and_recoveries(data, region_key = region_key, release_ndx_to_plot = 1:90, release_region_to_plt = "EGOA")
EGOA_releases
ggsave(filename = file.path(fig_path, "EGOA_tag_release_and_recovery.png"), width = 12, height = 8)

rec_dev_yrs_not_to_est = c(data$years[1:10], data$years[(length(data$years) - 1):length(data$years)])

## estimate 
# turn off tag-phi because we don't have tag-likelihood component
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi")
## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = NULL,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = "off", # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = F, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)
plot_SSB(mle_report)
plot_index_fit(mle_report)
gr_pars = as.vector(mle_obj$gr(mle_spatial$par))
names(gr_pars) = names(mle_spatial$par)

round(gr_pars, 4)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}
######################
## Include regional recruitment
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_NoTagData_regrecruit_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
data$global_rec_devs = 0
if(data$rec_devs_sum_to_zero == 1) {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
    
  }
} else {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = data$n_regions, ncol =  n_years)
  }
}

## check parameters
validate_input_data_and_parameters(data, parameters)

## save initial data and parameters
saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))

## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = NULL,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = "off", # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = F, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))

## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)
plot_SSB(mle_report)
gr_pars = as.vector(mle_obj$gr(mle_spatial$par))
names(gr_pars) = names(mle_spatial$par)

round(gr_pars, 4)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}
######################
## Estimate movement without tagging data
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_NoTagData_est_move_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

data$apply_fixed_movement = 0
## check parameters
validate_input_data_and_parameters(data, parameters)

## save initial data and parameters
saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))

## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = NULL,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = "off", # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = T, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))

## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_report_1 = mle_obj$report(mle_spatial$par)
plot_SSB(mle_report)
gr_pars = as.vector(mle_obj$gr(mle_spatial$par))
names(gr_pars) = names(mle_spatial$par)

round(gr_pars, 4)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Add movement matrix from Fenske's dissertation
######################
scenario = paste0("1960_fenske_movement_countrysrv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
data$apply_fixed_movement = 1
data$fixed_movement_matrix[,,1] = fenske_move_matrix / 100

## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = NULL,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = "off", # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = F, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))

## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

plot_SSB(mle_report, region_key = region_key)
plot_index_fit(mle_report, region_key = region_key)
tag_fits = get_tag_recovery_obs_fitted_values(mle_report, region_key = region_key)
tag_fits$residual = tag_fits$observed - tag_fits$predicted
tag_fits$recovery_region = factor(tag_fits$recovery_region, levels = region_key$area)
ggplot(tag_fits %>% filter(release_year %in% 1978:1985), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_1.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 1986:1995), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_2.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 1996:2005), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_3.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 2006:2015), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_4.png"), width = 10, height = 10)


sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Add movement matrix from Fenske's dissertation
## with spatail Q's
######################
scenario = paste0("1960_fenske_movement_spatial_q_countrysrv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
data$apply_fixed_movement = 1
data$fixed_movement_matrix[,,1] = fenske_move_matrix / 100

## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = NULL,  ## don't estimate the last one
                                   srv_q_spatial = T,
                                   tag_reporting_rate = "off", # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = F, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))

## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

plot_SSB(mle_report, region_key = region_key)
get_qs(mle_report, region_key = region_key, survey_labels = survey_labels)

plot_index_fit(mle_report, region_key = region_key)
tag_fits = get_tag_recovery_obs_fitted_values(mle_report, region_key = region_key)
tag_fits$residual = tag_fits$observed - tag_fits$predicted
tag_fits$recovery_region = factor(tag_fits$recovery_region, levels = region_key$area)
ggplot(tag_fits %>% filter(release_year %in% 1978:1985), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_1.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 1986:1995), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_2.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 1996:2005), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_3.png"), width = 10, height = 10)

ggplot(tag_fits %>% filter(release_year %in% 2006:2015), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(fig_path, "tag_fits_60_4.png"), width = 10, height = 10)


sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Add tag-recovery data in to the joint likelihood 
## assumed to be Poisson
## and estimate movement
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_TagData_01_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

# turn back on tag data
data$evaluate_tag_likelihood = 1
data$tag_likelihood = 0
# Estimate movement
data$apply_fixed_movement = 0

## some pars to fix

map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate, # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

plot_SSB(mle_report)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Add tag-recovery data in to the joint likelihood 
## assumed to be Poisson
## and estimate movement
## Fix second parameter for trawl selectivity between sexes
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_TagData_01_alt_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

read_in_previous_model_data_and_params = T
if(read_in_previous_model_data_and_params) {
  previous_scenario = paste0("1960_TagData_01_srv", ifelse(AF_direct_ageing, "a", ""))
  previous_fig_path = file.path("Output", "3-Area", previous_scenario)
  data = readRDS(file.path(previous_fig_path, "data.RDS"))
  parameters = readRDS(file.path(previous_fig_path, "parameters.RDS"))
  prev_mle = readRDS(file.path(previous_fig_path, "mle_report.RDS"))
}
# turn back on tag data
data$evaluate_tag_likelihood = 1
data$tag_likelihood = 0
# Estimate movement
data$apply_fixed_movement = 0
trwl_sel_second_shared_by_sex  = T
## some pars to fix

map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate, # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Change tag-recovery likelihood from 
## Poisson to Negative Binomial
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_TagData_02_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## Change likelihood switch
data$tag_likelihood = 1


## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate, # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Change tag-reporting rate to 
## be decadal
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_TagData_03_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## Decadal tag-reporting parameters, shared across regions
tag_report_time_block = c(1978, 1990, 2000, 2010)


## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   tag_reporting_rate = tag_report_time_block
                                   
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}

######################
## Change comp likelihood to Dirichlet-Multinomial 
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_Model_04_srv", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
## change likelihoods
data$fixed_catchatage_comp_likelihood = 1
data$trwl_catchatlgth_comp_likelihood = 1
data$fixed_catchatlgth_comp_likelihood = 0 ## had problems estimating over-dispersion parameter so setting to multinomial
parameters$trans_fixed_catchatlgth_error = log(1) # which was the estimated value for trawl LF
data$srv_catchatage_comp_likelihood = rep(1, data$n_surveys)

## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate, # turn off
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = F,
                                   est_trwl_LF_theta = T,
                                   est_srv_AF_theta = T
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

male_srv_AF = get_AF(MLE_report = start_report, observation = "srv", subset_years = 1979:2020, sex = "male", survey_labels = survey_labels)
ggplot(male_srv_AF %>% filter(Survey == "Japanese"), aes(x = Age)) +
  geom_point(aes(y = Observed, col = "Observed", shape = Sex, group = Sex)) +
  geom_line(aes(y = Predicted, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
  guides( linewidth = "none") +
  labs(y = "AF", col = "", linetype = "") +
  facet_grid(Year ~ Region) +
  theme_bw()
ggplot(male_srv_AF %>% filter(Survey == "Domestic"), aes(x = Age)) +
  geom_point(aes(y = Observed, col = "Observed", shape = Sex, group = Sex)) +
  geom_line(aes(y = Predicted, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
  guides( linewidth = "none") +
  labs(y = "AF", col = "", linetype = "") +
  facet_grid(Year ~ Region) +
  theme_bw()

plot_AF(MLE_report = start_report, observation = "srv", subset_years = 1979:1999, sex = "male")
plot_AF(start_report, observation = "srv", subset_years = 1979:1992, sex = "female")

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}


######################
## Spatial varying qs
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("1960_spatial_qs", ifelse(AF_direct_ageing, "a", ""))
fig_path = file.path("Output", "3-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
tag_report_time_block = c(1978, 1990, 2000, 2010)

read_in_previous_model_data_and_params = T
if(read_in_previous_model_data_and_params) {
  previous_scenario = paste0("1960_TagData_03_srv", ifelse(AF_direct_ageing, "a", ""))
  previous_fig_path = file.path("Output", "3-Area", previous_scenario)
  data = readRDS(file.path(previous_fig_path, "data.RDS"))
  parameters = readRDS(file.path(previous_fig_path, "parameters.RDS"))
  prev_mle = readRDS(file.path(previous_fig_path, "mle_report.RDS"))
}
## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = T,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement, # turn on and off to estimate
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   tag_reporting_rate = tag_report_time_block
                                   
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)
## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))


## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()

DEBUG_estimation = F
if(DEBUG_estimation) {
  mle_obj$env$tracepar = F
  options(warn=2) # covert warnings to errors and stop optimisations
}
## Optimisation
start_time = Sys.time()
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$estimation_time = Sys.time() - start_time
mle_spatial$evaluations
mle_spatial$convergence
mle_spatial$objective

## do an additional two Newton Raphson iterations to try and improve the fit.
## this will also check that the Hessian is well defined because it is needed
## to do the Newton Raphson iterations
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }
                       , error = function(e){e})

try_improve

## do some post-optimisation 
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)

## have a look at parameters and derived quantities
mle_report = mle_obj$report(mle_spatial$par)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_spatial$par)
saveRDS(mle_param_list, file.path(fig_path, "mle_par_list.RDS"))

## simulate observations
if(redo_sim_obs) {bl
  sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
  ## save this so we don't have to keep running it
  saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))
}
