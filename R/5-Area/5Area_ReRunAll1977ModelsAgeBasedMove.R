#'
#' Re-run all 5-Area model permutations
#' The other scripts are useful for doing a more detailed investigation on a model
#' This is to make sure that when they are all re-run they are
#' consistent in assumptions other than those that are specifically investigated
#' I found when I was doing this over multiple scripts that subtle things would change that will
#' cause inconsistencies between model runs other than those that are explicitly being investigated, which caused confusion
#' when summarising results
#'
#'
#' Notes: Changes Matt made: 
#' 1) Changed variable names to match the .cpp parameters (i.e., srv_dom_ll to srv)
#' 2) Added in some parameters that were not in here previously
#' 3) Changed the time block structure to be more consistent with other runs

# source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
library(here)
source(here("R", "Utility_Fxns.R"))
dat_path <- here("Data", "5-Area")

AF_direct_ageing = T

max_N_eff = 500 / 5 ## during the initial input calculations some of the input sample sizes can be quite large, which I don't want effecting the model
min_N_eff = 5 / 50
N_eff_multiplier = 0.4 ## only relevant for mult

## read in observational datasets
fixed_gear_AF_alk_pooled = readRDS(file = file.path(dat_path, "Observer_ALK_AF_w_eff.RDS"))
fixed_gear_AF_direct = readRDS(file = file.path(dat_path, "Observer_direct_AF_w_eff.RDS"))
trawl_observer_LF = readRDS(file = file.path(dat_path, "Observer_trawl_LF_w_eff.RDS"))
fixed_observer_LF = readRDS(file = file.path(dat_path, "Observer_fixed_LF_w_eff.RDS"))
survery_AF_direct = readRDS(file = file.path(dat_path, "Survey_direct_AF_w_eff.RDS")) %>% mutate(country = ifelse(year <= 1993, "Japan", "United States"))
survery_AF_alk_pooled = readRDS(file = file.path(dat_path, "Survey_ALK_AF_w_eff.RDS")) %>% mutate(country = ifelse(year <= 1993, "Japan", "United States"))

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

## survey index
design_survey_index = readRDS(file = file.path("Data", "Survey", "survey_ndx_5_area_combined_countries.RDS"))
## compare indicies between the three model structures for the Tweedie likelihood
design_survey_index$model = "Design"

## Tag data
tag_recovery_df = readRDS(file = file.path(dat_path, "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = file.path(dat_path, "Tag_release_summarised.RDS"))

## read in Catch
full_catch_df = readRDS(file = file.path(dat_path, "Catch_by_year_area_gear.RDS"))

## bring ADMB inputs to help configure initial model
sab_curr <- dget(file.path("Data", "ADMB", "tem.rdat")) 
sab_rep <-readLines(file.path("Data", "ADMB", "sable.rep"))
sab_inputs <- digest_sab_input(sablefish_input_filename = file.path("Data", "ADMB", "tem_2022_na_wh.dat")) 
sab_par <- readLines(file.path("Data", "ADMB", "tem.par")) 
sab_ctl <- digest_sab_ctl(sablefish_control_filename = file.path("Data", "ADMB", "tem.ctl")) 
## Model dimensions

years = 1977:2020
min_age = 2
max_age = 31
ages = min_age:max_age
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))
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
data$model = "TagIntegratedAgeBasedMovement"
data$ages = ages
data$years = years
data$length_bins = length_bins
data$n_regions = nrow(region_key)
survey_labels = "Alaska LL"
data$n_surveys = length(survey_labels)
n_regions = data$n_regions
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
data$n_init_rec_devs = 15 ## number of non-equilibrium age-structure
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
move_matrix = data$fixed_movement_matrix_old = data$fixed_movement_matrix_young = matrix(0, nrow = n_regions, ncol = n_regions);
diag(move_matrix) = diag(data$fixed_movement_matrix_old) = diag(data$fixed_movement_matrix_young) = 1

move_matrix = move_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1)
# renormalise
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/")

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
data$prop_F_hist = 1
# drop years outside model years
full_catch_df = full_catch_df %>% filter(year %in% years)
full_catch_df$year_f = factor(full_catch_df$year, levels = year_lvls, ordered = T)
full_catch_df$region_f = factor(full_catch_df$area, levels = reg_lvls, ordered = T)
## convert catch to metric tonnes
full_catch_df$catch_mt = full_catch_df$catch / 1000
trwl_catch = full_catch_df %>% filter(fmp_gear == "TRW") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)
fixed_catch = full_catch_df %>% filter(fmp_gear == "HAL") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)

## replace 0 values with low values
fixed_catch = fixed_catch %>% dplyr::select(-region_f) %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))
trwl_catch = trwl_catch %>% dplyr::select(-region_f) %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))

data$trwl_fishery_catch = as.matrix(trwl_catch, nrow = n_regions, ncol = n_projyears)
data$fixed_fishery_catch = as.matrix(fixed_catch, nrow = n_regions, ncol = n_projyears)

## control variables
# 3 time-blocks for fixed gear selectivity
data$fixed_sel_type = as.vector(rep(0, 3), mode = "integer")
data$fixed_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")
# data$fixed_sel_by_year_indicator[projyears %in% 1995:2015] = 1
# data$fixed_sel_by_year_indicator[projyears > 2015] = 2
data$fixed_sel_by_year_indicator[projyears > 2015] = 1


# single time-block for trawl fishery
data$trwl_sel_type = as.vector(rep(1, 1), mode = "integer")
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

## two time-blocks for survey selectivity (not sure why... )
data$srv_dom_ll_sel_type = as.vector(rep(0, 2), mode = "integer")
data$srv_dom_ll_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")
# data$srv_dom_ll_sel_by_year_indicator[projyears %in% 1988:1994] = 1
# data$srv_dom_ll_sel_by_year_indicator[projyears > 1994] = 2

data$srv_sel_type = matrix(0, nrow = 1, ncol = data$n_surveys)
data$srv_sel_by_year_indicator = matrix(0, nrow = n_projyears, ncol = data$n_surveys)
# data$srv_sel_by_year_indicator[projyears %in% 1988:1994,] = 1
# data$srv_sel_by_year_indicator[projyears > 1994,] = 2

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

data$n_years_to_retain_tagged_cohorts_for = 10
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
  obs_years = obs_years[!obs_years %in% unique(fixed_gear_AF$year)]
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
data$srv_dom_ll_catchatage_indicator = as.matrix(survey_AF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$srv_catchatage_indicator = array(as.matrix(survey_AF_indicator %>% dplyr::select(-region_f)), dim = c(n_regions, n_years,data$n_surveys))
data$obs_srv_dom_ll_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years))
data$obs_srv_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years, data$n_surveys))
data$srv_bio_likelihood = rep(1, data$n_surveys)
data$srv_q_transformation = rep(0, data$n_surveys)

obs_years = unique(survery_AF$year)
obs_reg = unique(survery_AF$region)
## populate container
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (survery_AF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx])) %>% arrange(sex_age)
    if(nrow(this_df) > 0) {
      data$obs_srv_catchatage[,this_region_ndx, this_year_ndx,1] = this_df$P * this_df$eff_N
      data$obs_srv_dom_ll_catchatage[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N
    }
  }
}
data$srv_dom_ll_catchatage_covar_structure = 0
data$srv_dom_ll_catchatage_comp_likelihood = 0
data$srv_catchatage_covar_structure = 0
data$srv_catchatage_comp_likelihood = 0

### Survey LL index
data$srv_dom_ll_bio_indicator = as.matrix(as.matrix(survey_index_indicator %>% dplyr::select(-region_f)), nrow = n_regions, ncol = n_years)
data$srv_bio_indicator = array(as.matrix(survey_index_indicator %>% dplyr::select(-region_f)), dim = c(n_regions, n_years, data$n_surveys))
data$obs_srv_dom_ll_bio = array(0, dim = c(n_regions, n_years))
data$obs_srv_bio = array(0, dim = c(n_regions, n_years, data$n_surveys ))
data$obs_srv_dom_ll_se = array(0.0, dim = c(n_regions, n_years))
data$obs_srv_se = array(0, dim = c(n_regions, n_years, data$n_surveys ))


obs_years = unique(survey_index$Year)
obs_reg = unique(survey_index$area_lab)
likelihood_type = "Lognormal"
## populate container
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (survey_index %>% filter(area_lab == obs_reg[r_ndx], Year == obs_years[y_ndx]))
    if(nrow(this_df) > 0) {
      data$obs_srv_dom_ll_bio[this_region_ndx, this_year_ndx] = this_df$sum_estimates## model numbers are in 1000's is in kilo tonnes
      data$obs_srv_bio[this_region_ndx, this_year_ndx,1] = this_df$sum_estimates## model numbers are in 1000's is in kilo tonnes
      
      #data$obs_srv_dom_ll_se[this_region_ndx, this_year_ndx] = this_df$se
      data$obs_srv_dom_ll_se[this_region_ndx, this_year_ndx] = (this_df$se)  # scale by 1000 as well
      data$obs_srv_se[this_region_ndx, this_year_ndx,1] = (this_df$se)  # scale by 1000 as well
      
    }
  }
}

data$srv_dom_ll_bio_likelihood = 0
data$srv_dom_ll_obs_is_abundance = 1
data$srv_obs_is_abundance = rep(1, data$n_surveys)
data$srv_dom_ll_q_by_year_indicator = matrix(0, nrow = n_years, ncol = data$n_surveys)
## three q time-blocks
# data$srv_dom_ll_q_by_year_indicator[projyears %in% 1988:1994,] = 1
# data$srv_dom_ll_q_by_year_indicator[projyears > 1994,] = 2

data$srv_q_by_year_indicator = matrix(0, nrow = n_years, ncol = data$n_surveys)

## three q time-blocks? (not sure why since the other models dont really do that)
# data$srv_q_by_year_indicator[projyears %in% 1988:1994,] = 1
# data$srv_q_by_year_indicator[projyears > 1994,] = 2

data$srv_dom_ll_q_transformation = 0
data$q_is_nuisance = 0

###
# Tag-recoveries
# This can take a bit of time to get your head around!
###
include_tag_recoveries = T
include_zero_tag_recovery_events = T
tag_recovery_years = 1978:2016

data$tag_recovery_indicator_by_year = rep(0, n_years) ## no tag releases
data$obs_tag_recovery = array(0, dim = c(n_ages, n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))
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
      recovery_region_ndx = region_key$TMB_ndx[which(region_key$area %in% recovery_regions[r_ndx])] + 1
      ## tags must be at liberty for at least one year
      this_recovery_year_region = this_recovery_year %>% filter(release_year <= tag_recovery_years[y_ndx] - 1,release_year %in% tag_release_years, recovery_region == recovery_regions[r_ndx])
      # now find how many 'release-event' tags we recovered during this recovery year and region.- which is release year x release region specific blahh!!!
      for(release_event_ndx in 2:(data$n_years_to_retain_tagged_cohorts_for)) { # index starts at two because 1 indicates they were released this year and we don't consider recoveries until after a year at liberty
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
          release_ndx = region_key$TMB_ndx[which(region_key$area %in% release_regions[rel_ndx])] + 1
          
          mod_release_ndx = which(tag_release_years %in% this_release_year)
          if((sum(data$male_tagged_cohorts_by_age[,release_ndx, mod_release_ndx]) + sum(data$female_tagged_cohorts_by_age[,release_ndx, mod_release_ndx])) <= 0)
            next;
          model_tag_release_ndx = get_tag_release_ndx(release_ndx, release_event_ndx, n_regions)
          ## there was a possible observation for this release event
          if(include_zero_tag_recovery_events)
            data$tag_recovery_indicator[model_tag_release_ndx, recovery_region_ndx, y_ndx] = 1
          ## not all regions were 'release events' check we actually released fish in this event
          specific_recovery_df = this_release_event_df %>% filter(region_release == release_regions[rel_ndx])
          if(nrow(specific_recovery_df) > 0) {
            ## tell the model a recovery observation happened for this release event
            data$tag_recovery_indicator[model_tag_release_ndx, recovery_region_ndx, y_ndx] = 1
            ## get the recovery age frequency
            data$obs_tag_recovery[,model_tag_release_ndx, recovery_region_ndx, y_ndx] = ((specific_recovery_df %>% group_by(age) %>% summarise(Nage_at_recovery = sum(Nage_at_recovery)))$Nage_at_recovery)
          }
        }
      }
    }
  }
}
## are there tag-recoveries in the first year at liberty. There shouldn't be
if(sum(data$tag_recovery_indicator[1:data$n_regions,,]) != 0)
  cat("found recoveries in the first year of release\n")
if(sum(data$tag_recovery_indicator[(dim(data$tag_recovery_indicator)[1] - data$n_regions + 1):dim(data$tag_recovery_indicator)[1],,]) != 0)
  cat("found recoveries in the plus group\n")

cat("The number of recovery observations ", sum(data$tag_recovery_indicator), "\n")
data$tag_likelihood = 2
data$evaluate_tag_likelihood = 1
data$age_based_movement = 1
data$movement_time_block_indicator = rep(0, n_years)
data$n_movement_time_blocks = 1


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
# parameters$ln_srv_dom_ll_sel_pars = array(0, dim = c(length(unique(data$srv_dom_ll_sel_by_year_indicator )), 2, 2))
parameters$ln_srv_sel_pars = array(0, dim = c(length(unique(data$srv_dom_ll_sel_by_year_indicator )), 2, 2))

## populate parameters Note some of the male delta values are set to the female values. Line 1800 tem.tpl

## populate parameters Note some of the male delta values are set to the female values. Line 1800 tem.tpl
if(dim(parameters$ln_fixed_sel_pars)[1] == 2) {
  parameters$ln_fixed_sel_pars[1,1,1]  = parameters$ln_fixed_sel_pars[2,1,1]  =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  
} else {
  parameters$ln_fixed_sel_pars[1,1,1]  = parameters$ln_fixed_sel_pars[2,1,1]  = parameters$ln_fixed_sel_pars[3,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = parameters$ln_fixed_sel_pars[3,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = parameters$ln_fixed_sel_pars[3,1,2] =suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
  parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] = parameters$ln_fixed_sel_pars[3,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
}
parameters$ln_trwl_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_m:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))


## NOTE: all delta parameters for all survey selectivities are fixed based on srv_dom_ll_1
if(dim(parameters$ln_srv_sel_pars)[1] == 1) {
  parameters$ln_srv_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_sel_pars)[1] == 2) {
  parameters$ln_srv_sel_pars[1,1,1] =  parameters$ln_srv_sel_pars[2,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1] = parameters$ln_srv_sel_pars[2,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2] = parameters$ln_srv_sel_pars[2,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2] = parameters$ln_srv_sel_pars[2,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_sel_pars)[1] == 3) {
  parameters$ln_srv_sel_pars[1,1,1] =parameters$ln_srv_sel_pars[2,1,1] = parameters$ln_srv_sel_pars[3,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1] = parameters$ln_srv_sel_pars[2,2,1] = parameters$ln_srv_sel_pars[3,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2] = parameters$ln_srv_sel_pars[2,1,2] = parameters$ln_srv_sel_pars[3,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2] = parameters$ln_srv_sel_pars[2,2,2] = parameters$ln_srv_sel_pars[3,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
}

parameters$ln_srv_sel_pars = array(parameters$ln_srv_sel_pars, dim = c(length(unique(data$srv_dom_ll_sel_by_year_indicator )), 2, 2, 1))


## movement pars
parameters$transformed_movement_pars_young = matrix(NA, nrow = n_regions - 1, ncol = n_regions)
for(i in 1:n_regions)
  parameters$transformed_movement_pars_young[,i] = simplex(move_matrix[i,])
parameters$transformed_movement_pars_old = parameters$transformed_movement_pars_young

# set to maturity initially
plot(data$ages, data$maturity[,1], type = "l", lwd = 3, xlab = "ages", ylab = "Proportion mature")
lines(data$ages, logis(data$ages, 7.05, 4.0), lty = 2, col = "red", lwd = 3)
parameters$ln_ato95_movement = log(4.0)
parameters$ln_a50_movement = log(7.05)

parameters$ln_fixed_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish1:",sab_par)+1],split=" "))))
parameters$ln_fixed_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish1:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_trwl_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish3:",sab_par)+1],split=" "))))
parameters$ln_trwl_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish3:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_init_F_avg = parameters$ln_fixed_F_avg
# parameters$trans_srv_dom_ll_q = array(log(7), dim =c(n_regions, length(unique(data$srv_dom_ll_q_by_year_indicator))))
parameters$trans_srv_q = array(log(7), dim =c(n_regions, length(unique(data$srv_dom_ll_q_by_year_indicator)), data$n_surveys))

admb_ln_rec_devs = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_rec_dev:",sab_par)+1],split=" "))))
admb_ln_rec_devs = subset(admb_ln_rec_devs, subset = !is.na(admb_ln_rec_devs))


if(data$rec_devs_sum_to_zero == 1) {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
    
  }
} else {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years)
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
parameters$logistic_prop_recruit_male = rep(logit(0.5), length(data$years))
parameters$trans_SR_pars = rep(qlogis(0.8), 1)
parameters$trans_srv_catchatage_error = rep(log(1), data$n_surveys)



#'
#' Initial parameter setup that will be shared among 
#' the model runs
#' these relate to parameters that are passed to the
#' set_up_parameters function, which controls what parameters are estimated
#' or shared. This was done to ensure these were consistent among multiple runs
#' unless a run was specifically investigating/changing this assumption

srv_sel_first_param_shared_by_sex = T
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_shared_by_sex  = T
fixed_sel_second_shared_by_sex   = T
trwl_sel_first_shared_by_sex  = T
trwl_sel_second_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021  ## don't estimate the last one
srv_q_spatial = F
tag_reporting_rate = "constant"
est_init_F = F
est_catch_sd = F
est_movement = T
est_prop_male_recruit = "off"


######################
## First model run 01a
######################
# options(warn=0) ## if warnings are being set to errors
scenario = paste0("Model1977_01", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

data$age_based_movement = 0
data$tag_likelihood = 0 # mutlinomial...
## check parameters
validate_input_data_and_parameters(data, parameters)

## save initial data and parameters
saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))

## plot some input info
plot_input_observations(data, region_key = region_key)
ggsave(filename = file.path(fig_path, "Observation_Frequency.png"), width = 10, height = 10)

plot_input_timeblocks(data)
ggsave(filename = file.path(fig_path, "time_blocks.png"), width = 7, height = 6)


plot_input_catches(data, region_key = region_key)
ggsave(filename = file.path(fig_path, "InputCatches.png"), width = 8, height = 6)

## estimate 
na_map = fix_pars(par_list = parameters,
                  pars_to_exclude = c("ln_ato95_movement", "ln_a50_movement", "ln_tag_phi"))
## some pars to fix
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = na_map,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data = data, parameters = parameters, map = map_fixed_pars, 
                     DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = F)

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
mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, 
                     control = list(iter.max = 10000, eval.max = 10000))
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
plot(data$ages, mle_report$old_age_based_movement_ogive, type = "l", lwd = 3, col = "black")
lines(data$ages, mle_report$old_age_based_movement_ogive, lty = 2, lwd = 3, col = "red")

mle_report$a50_movement
mle_report$ato95_movement

round(mle_report$movement_matrix_old, 2)
round(mle_report$movement_matrix_young, 2)

sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))


######################
## Run model run 02a
## Change composition likelihoods to
## Dirichlet-multinomial
######################
scenario = paste0("Model1977_02", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## change likelihoods
data$fixed_catchatage_comp_likelihood = 1
data$trwl_catchatlgth_comp_likelihood = 1
data$fixed_catchatlgth_comp_likelihood = 1
data$srv_dom_ll_catchatage_comp_likelihood = 1

## check they are consistent
validate_input_data_and_parameters(data, parameters)

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
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))


######################
## Run model run 03a
## Change recruitment deviations to be 
## regional instead of global
######################
scenario = paste0("Model1977_03", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## Change recruitment deviation to be regional
data$global_rec_devs = 0
if(data$rec_devs_sum_to_zero == 1) {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
    
  }
} else {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years)
  }
}

## check they are consistent
validate_input_data_and_parameters(data, parameters)

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
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
                                   est_trwl_LF_theta = T,
                                   est_srv_ll_AF_theta = T
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))


######################
## Run model run 04a
## force regional recruitment deviations to sum = 0
######################
scenario = paste0("Model1977_04", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## Change regional recruitment deviations to sum = 0
data$rec_devs_sum_to_zero = 1
if(data$rec_devs_sum_to_zero == 1) {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
    
  }
} else {
  if(data$global_rec_devs == 1) {
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years)
  }
}

## check they are consistent
validate_input_data_and_parameters(data, parameters)

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
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
                                   est_trwl_LF_theta = T,
                                   est_srv_ll_AF_theta = T
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))

######################
## Run model run 05a
## Change tag-likelihood to be Poisson instead of Negative Binomial
######################
scenario = paste0("Model1977_05", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

# change tag-likelihood
data$tag_likelihood = 0

## check they are consistent
validate_input_data_and_parameters(data, parameters)

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
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
                                   est_trwl_LF_theta = T,
                                   est_srv_ll_AF_theta = T
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))


######################
## Run model run 06a
## Estimate tag-reporting rates for every decade
######################
scenario = paste0("Model1977_06", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

tag_report_time_block = c(1978, 1990, 2000, 2010)
## this is used in the function 'set_up_parameters'

## check they are consistent
validate_input_data_and_parameters(data, parameters)

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
                                   tag_reporting_rate = tag_report_time_block,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
                                   est_trwl_LF_theta = T,
                                   est_srv_ll_AF_theta = T
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))

######################
## Run model run 07
## Estimate decadal tag-reporting rates 
## with negative binomial
## Change trawl selectivity to three parameter double normal
######################
scenario = paste0("Model1977_07", ifelse(AF_direct_ageing, "a", ""), "_AgeBasedMove")
fig_path = file.path("Output", "5-Area", scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

read_in_previous_model_data_and_params = T
if(read_in_previous_model_data_and_params) {
  previous_scenario = paste0("Model1977_06", ifelse(AF_direct_ageing, "a", ""))
  previous_fig_path = file.path("Output", "5-Area", previous_scenario)
  data = readRDS(file.path(previous_fig_path, "data.RDS"))
  parameters = readRDS(file.path(previous_fig_path, "parameters.RDS"))
  
}

tag_report_time_block = c(1978, 1990, 2000, 2010)
## Change likelihood to Negative Binomial
data$tag_likelihood = 1
# Change selectivity to three parameter double normal
data$trwl_sel_type = as.vector(rep(5, 1), mode = "integer")
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

parameters$ln_trwl_sel_pars = array(log(5), dim = c(1,3,2))
parameters$ln_trwl_sel_pars[1,1,] = log(8) # mean
parameters$ln_trwl_sel_pars[1,2,] = log(5) # sig r
parameters$ln_trwl_sel_pars[1,3,] = log(2) # sig l


## check they are consistent
validate_input_data_and_parameters(data, parameters)

map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = F,
                                   trwl_sel_second_shared_by_sex  = T, # estimate a common sigma r between sexes
                                   trwl_sel_third_shared_by_sex  = T, # estimate a common sigma l between sexes
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,  ## don't estimate the last one
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_report_time_block,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit,
                                   est_fixed_AF_theta = T,
                                   est_fixed_LF_theta = T,
                                   est_trwl_LF_theta = T,
                                   est_srv_ll_AF_theta = T
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

## simulate observations
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))



