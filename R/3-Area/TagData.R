#'
#' The first attempt at fitting data to the spatial model
#' with tagging in the partition
#'

source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
# options(warn=0) ## if warnings are being set to errors
scenario = "TagData_01"
fig_path = file.path(DIR$app_3A, scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}


## read in observational datasets
fixed_gear_AF_alk_pooled = readRDS(file = file.path(DIR$data3A, "Observer_ALK_AF_w_eff.RDS"))
fixed_gear_AF_direct = readRDS(file = file.path(DIR$data3A, "Observer_direct_AF_w_eff.RDS"))
trawl_observer_LF = readRDS(file = file.path(DIR$data3A, "Observer_trawl_LF_w_eff.RDS"))
fixed_observer_LF = readRDS(file = file.path(DIR$data3A, "Observer_fixed_LF_w_eff.RDS"))
survery_AF_direct = readRDS(file = file.path(DIR$data3A, "Survey_direct_AF_w_eff.RDS"))
survery_AF_alk_pooled = readRDS(file = file.path(DIR$data3A, "Survey_ALK_AF_w_eff.RDS"))
max_N_eff = 1000 ## during the initial input calculations some of the input sample sizes were huge, which I don't want effecting the model
min_N_eff = 50

N_eff_multiplier = 0.4

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

design_survey_index = readRDS(file = file.path(DIR$data, "Survey", "survey_ndx_5_area_combined_countries.RDS"))
design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "AI" ~ "BS_AI_WGOA",
                                                                 area_lab == "BS" ~ "BS_AI_WGOA",
                                                                 area_lab == "WGOA" ~ "BS_AI_WGOA",
                                                                 TRUE ~ area_lab))
## sum estimates for BS, AI, WGOA
design_survey_index = design_survey_index %>% group_by(Year, area_lab) %>% summarise(sum_estimates = sum(sum_estimates), sum_var = sum(sum_var), se = sqrt(sum_var), LCI = sum_estimates - 2*se, UCI = sum_estimates + 2*se)

ggplot(design_survey_index, aes(Year, sum_estimates, col = area_lab, linetype = area_lab, fill = area_lab, alpha = 0.3)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI)) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Relative index (RPN)", linetype = "", col = "", fill = "") +
  guides(alpha = "none") +
  facet_wrap(~area_lab, scales = "free_y") +
  theme_bw()
ggsave(filename = file.path(DIR$input_figs_3A, "indicies.png"), width = 7, height = 7)

## Tag data
tag_recovery_df = readRDS(file = file.path(DIR$data3A, "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = file.path(DIR$data3A, "Tag_release_summarised.RDS"))

## read in Catch
full_catch_df = readRDS(file = file.path(DIR$data3A, "Catch_by_year_area_gear.RDS"))

## bring ADMB inputs to help configure initial model
sab_curr <- dget(file.path(DIR$admb, "tem.rdat")) 
sab_rep <-readLines(file.path(DIR$admb, "sable.rep"))
sab_inputs <- digest_sab_input(sablefish_input_filename = file.path(DIR$admb, "tem_2022_na_wh.dat")) 
sab_par <- readLines(file.path(DIR$admb, "tem.par")) 
sab_ctl <- digest_sab_ctl(sablefish_control_filename = file.path(DIR$admb, "tem.ctl")) 
## Model dimensions

years = 1977:2020
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
data$ages = ages
data$years = years
data$length_bins = length_bins
data$n_regions = nrow(region_key)
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
data$n_init_rec_devs = 28 ## number of non-equilibrium age-structure
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
data$movement_matrix = data$fixed_movement_matrix = matrix(0, nrow = n_regions, ncol = n_regions);
diag(data$movement_matrix) = 0.9
diag(data$fixed_movement_matrix) = 1
data$movement_matrix = data$movement_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1)
# renormalise
data$movement_matrix = sweep(data$movement_matrix, 1, STATS = rowSums(data$movement_matrix), "/")

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
data$fixed_sel_by_year_indicator[projyears %in% 1995:2015] = 1
data$fixed_sel_by_year_indicator[projyears > 2015] = 2

# single time-block for trawl fishery
data$trwl_sel_type = as.vector(rep(1, 1), mode = "integer")
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

## two time-blocks for survey selectivity
data$srv_dom_ll_sel_type = as.vector(rep(0, 3), mode = "integer")
data$srv_dom_ll_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")
data$srv_dom_ll_sel_by_year_indicator[projyears %in% 1988:1994] = 1
data$srv_dom_ll_sel_by_year_indicator[projyears > 1994] = 2

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
fixed_gear_AF_alk_pooled = fixed_gear_AF_alk_pooled %>% filter(year %in% years)
fixed_gear_AF_alk_pooled$region_f = factor(fixed_gear_AF_alk_pooled$region, levels = reg_lvls, ordered = T)
fixed_gear_AF_alk_pooled$year_f = factor(fixed_gear_AF_alk_pooled$year, levels = year_lvls, ordered = T)
trawl_observer_LF = trawl_observer_LF %>% filter(year %in% years)
trawl_observer_LF$region_f = factor(trawl_observer_LF$region, levels = reg_lvls, ordered = T)
trawl_observer_LF$year_f = factor(trawl_observer_LF$year, levels = year_lvls, ordered = T)
fixed_observer_LF = fixed_observer_LF %>% filter(year %in% years)
fixed_observer_LF$region_f = factor(fixed_observer_LF$region, levels = reg_lvls, ordered = T)
fixed_observer_LF$year_f = factor(fixed_observer_LF$year, levels = year_lvls, ordered = T)

survery_AF_alk_pooled = survery_AF_alk_pooled %>% filter(year %in% years)
survery_AF_alk_pooled$region_f = factor(survery_AF_alk_pooled$region, levels = reg_lvls, ordered = T)
survery_AF_alk_pooled$year_f = factor(survery_AF_alk_pooled$year, levels = year_lvls, ordered = T)
## drop years outside of model years. can cause NA's
survey_index = design_survey_index %>% filter(Year %in% years)
survey_index$region_f = factor(survey_index$area_lab, levels = reg_lvls, ordered = T)
survey_index$year_f = factor(survey_index$Year, levels = year_lvls, ordered = T)

indicator_fun = function(x) {
  ifelse(is.na(sum(x)), 0, 1)
}

fixed_AF_indicator = fixed_gear_AF_alk_pooled %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
trawl_LF_indicator = trawl_observer_LF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
fixed_LF_indicator = fixed_observer_LF %>% filter(!year %in% unique(fixed_gear_AF_alk_pooled$year)) %>% ungroup() %>% dplyr::select(region_f, year_f, P) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(P), values_fn = indicator_fun)
survey_AF_indicator = survery_AF_alk_pooled %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
survey_index_indicator = survey_index %>% ungroup() %>% dplyr::select(region_f, year_f, Year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(Year), values_fn = indicator_fun)

### Fixed gear fishery AF
data$fixed_catchatage_indicator = as.matrix(fixed_AF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$obs_fixed_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years))
## fill the container
obs_years = unique(fixed_gear_AF_alk_pooled$year)
obs_reg = unique(fixed_gear_AF_alk_pooled$region)
N_eff = 100 ## sample size
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (fixed_gear_AF_alk_pooled %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx]))
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
data$srv_dom_ll_catchatage_indicator = as.matrix(survey_AF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$obs_srv_dom_ll_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years))
obs_years = unique(survery_AF_alk_pooled$year)
obs_reg = unique(survery_AF_alk_pooled$region)
## populate container
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (survery_AF_alk_pooled %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx])) %>% arrange(sex_age)
    if(nrow(this_df) > 0)
      data$obs_srv_dom_ll_catchatage[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N
  }
}
data$srv_dom_ll_catchatage_covar_structure = 0
data$srv_dom_ll_catchatage_comp_likelihood = 0

### Survey LL index
data$srv_dom_ll_bio_indicator = as.matrix(survey_index_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years)
data$obs_srv_dom_ll_bio = array(0, dim = c(n_regions, n_years))
data$obs_srv_dom_ll_se = array(0.0, dim = c(n_regions, n_years))

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
      data$obs_srv_dom_ll_bio[this_region_ndx, this_year_ndx] = this_df$sum_estimates  ## model numbers are in 1000's is in kilo tonnes
      
      #data$obs_srv_dom_ll_se[this_region_ndx, this_year_ndx] = this_df$se
      data$obs_srv_dom_ll_se[this_region_ndx, this_year_ndx] = (this_df$se)  # scale by 1000 as well
    }
  }
}

data$srv_dom_ll_bio_likelihood = 0
data$srv_dom_ll_obs_is_abundance = 0
data$srv_dom_ll_q_by_year_indicator = rep(0, n_years)
## three q time-blocks
data$srv_dom_ll_q_by_year_indicator[projyears %in% 1988:1994] = 1
data$srv_dom_ll_q_by_year_indicator[projyears > 1994] = 2
data$srv_dom_ll_q_transformation = 0
data$q_is_nuisance = 0

###
# Tag-recoveries
# This can take a bit of time to get your head around!
###
include_tag_recoveries = T
include_zero_tag_recovery_events = T
tag_recovery_years = 1979:2020

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
data$evaluate_tag_likelihood = 1

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
parameters$ln_srv_dom_ll_sel_pars = array(0, dim = c(length(unique(data$srv_dom_ll_sel_by_year_indicator )), 2, 2))

## populate parameters Note some of the male delta values are set to the female values. Line 1800 tem.tpl

parameters$ln_fixed_sel_pars[1,1,1]  = parameters$ln_fixed_sel_pars[2,1,1]  = parameters$ln_fixed_sel_pars[3,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = parameters$ln_fixed_sel_pars[3,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = parameters$ln_fixed_sel_pars[3,1,2] =suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] = parameters$ln_fixed_sel_pars[3,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
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

## NOTE: all delta parameters for all survey selectivities are fixed based on srv_dom_ll_1
if(dim(parameters$ln_srv_dom_ll_sel_pars)[1] == 1) {
  parameters$ln_srv_dom_ll_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,1] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,1,2] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_dom_ll_sel_pars)[1] == 2) {
  parameters$ln_srv_dom_ll_sel_pars[1,1,1] =  parameters$ln_srv_dom_ll_sel_pars[2,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,1] = parameters$ln_srv_dom_ll_sel_pars[2,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,1,2] = parameters$ln_srv_dom_ll_sel_pars[2,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,2] = parameters$ln_srv_dom_ll_sel_pars[2,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
} else if(dim(parameters$ln_srv_dom_ll_sel_pars)[1] == 3) {
  parameters$ln_srv_dom_ll_sel_pars[1,1,1] =parameters$ln_srv_dom_ll_sel_pars[2,1,1] = parameters$ln_srv_dom_ll_sel_pars[3,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,1] = parameters$ln_srv_dom_ll_sel_pars[2,2,1] = parameters$ln_srv_dom_ll_sel_pars[3,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,1,2] = parameters$ln_srv_dom_ll_sel_pars[2,1,2] = parameters$ln_srv_dom_ll_sel_pars[3,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_dom_ll_sel_pars[1,2,2] = parameters$ln_srv_dom_ll_sel_pars[2,2,2] = parameters$ln_srv_dom_ll_sel_pars[3,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
}

## movement pars
parameters$transformed_movement_pars = matrix(NA, nrow = n_regions - 1, ncol = n_regions)
for(i in 1:n_regions)
  parameters$transformed_movement_pars[,i] = simplex(data$movement_matrix[i,])

parameters$ln_fixed_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish1:",sab_par)+1],split=" "))))
parameters$ln_fixed_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish1:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_trwl_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish3:",sab_par)+1],split=" "))))
parameters$ln_trwl_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish3:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5))))

parameters$ln_init_F_avg = parameters$ln_fixed_F_avg
parameters$trans_srv_dom_ll_q = array(log(7), dim =c(n_regions, length(unique(data$srv_dom_ll_q_by_year_indicator))))
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
parameters$trans_srv_dom_ll_catchatage_error = log(1)

parameters$logistic_prop_recruit_male = rep(logit(0.5), length(data$years))
######
## Check inputs
######
data$model = "TagIntegrated"

validate_input_data_and_parameters(data, parameters)

## save objects to speed up time in debugging
#
if(FALSE) {
  data = readRDS(file.path(fig_path, "data.RDS"))
  parameters = readRDS(file.path(fig_path, "parameters.RDS"))
}

## Create AD object
validate_input_data_and_parameters(data, parameters)

## plot some input info
plot_input_observations(data, region_key = region_key)
ggsave(filename = file.path(fig_path, "Observation_Frequency.png"), width = 10, height = 10)

plot_input_timeblocks(data)
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


#####################
## try and estimate 
#####################
## some pars to fix
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi")

map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = F,
                                   srv_sel_second_param_shared_by_sex = T,
                                   fixed_sel_first_shared_by_sex  = F,
                                   fixed_sel_second_shared_by_sex   = F,
                                   trwl_sel_first_shared_by_sex  = F,
                                   trwl_sel_second_shared_by_sex  = F,
                                   recruit_dev_years_not_to_estimate = NULL, 
                                   srv_q_spatial = F,
                                   tag_reporting_rate = "constant",
                                   est_init_F = T,
                                   est_catch_sd = F,
                                   est_movement = T,
                                   est_prop_male_recruit = "off"
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

data$apply_fixed_movement = 0
## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

## what parameters are being estimated
unique(names(mle_obj$par))

## pre-optim sanity checks
pre_optim_sanity_checks(mle_obj)
start_report = mle_obj$report()


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
## get standard errors
sd_report = sdreport(mle_obj)

saveRDS(data, file.path(fig_path, "data.RDS"))
saveRDS(parameters, file.path(fig_path, "parameters.RDS"))
saveRDS(mle_report, file.path(fig_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(fig_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(fig_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(fig_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(fig_path, "region_key.RDS"))

if(FALSE){
  data = readRDS(file.path(fig_path, "data.RDS"))
  region_key = readRDS(file.path(fig_path, "region_key.RDS"))
  parameters = readRDS(file.path(fig_path, "parameters.RDS"))
  mle_report = readRDS(file.path(fig_path, "mle_report.RDS"))
  sd_report = readRDS(file.path(fig_path, "sd_report.RDS"))
  mle_spatial = readRDS(file.path(fig_path, "mle_optim.RDS"))
  map_fixed_pars = readRDS(file.path(fig_path, "map_fixed_pars.RDS"))
}


nll_df = get_negloglike(mle_report)
sum(nll_df$negloglike)
## initial age-structure
plot_init_nage(mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "Init_numbers_at_age.png"), width = 8, height = 6)

## movement assumption
plot_movement(mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "Est_movement.png"), width = 8, height = 8)

# starting values
plot_movement(start_report, region_key)
# plot selectivities
plot_selectivities(mle_report)
ggsave(filename = file.path(fig_path, "Est_selectivity.png"), width = 6, height = 6)

# starting values
plot_selectivities(start_report)

## SSBs
ssb_plt = plot_SSB(mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "regional_ssb.png"), width = 8, height = 6)
plot_SSB(start_report, region_key = region_key)
ssb_plt

## sum SSB over all regions
ggplot(ssb_plt$data %>% group_by(Year) %>% summarise(total_ssb = sum(SSB)), aes(x = Year, y = total_ssb)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA)
ggsave(filename = file.path(fig_path, "global_ssb.png"), width = 8, height = 6)

## plot F's
plot_fishing_mortalities(MLE_report = mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "recuitment.png"), width = 8, height = 6)

## plot Recruitment
plot_recruitment(MLE_report = mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "recuitment.png"), width = 8, height = 6)

## plot partition
#plot_partition(MLE_report = mle_report, region_key = region_key, subset_years = 2015:2021)

####
# plot fits
####
## catch and index
plot_catch_fit(MLE_report = mle_report, region_key = region_key) + facet_wrap(label~Region, ncol = n_regions) + ylab("Catch (mt)")
ggsave(filename = file.path(fig_path, "Catch_fit.png"), width = 8, height = 7)
plot_index_fit(MLE_report = mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "survey_index_fit.png"), width = 8, height = 8)

## survey mean AFs
plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "srv_dom_ll",sex = "male")
ggsave(filename = file.path(fig_path, "survey_male_mean_AF.png"), width = 8, height = 6)
plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "srv_dom_ll",sex = "female")
ggsave(filename = file.path(fig_path, "survey_female_mean_AF.png"), width = 8, height = 6)

plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male")
ggsave(filename = file.path(fig_path, "fixed_fishery_male_mean_AF.png"), width = 8, height = 6)
plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female")
ggsave(filename = file.path(fig_path, "fixed_fishery_femlae_mean_AF.png"), width = 8, height = 6)

plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male")
ggsave(filename = file.path(fig_path, "fixed_fishery_male_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female")
ggsave(filename = file.path(fig_path, "fixed_fishery_femlae_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "male")
ggsave(filename = file.path(fig_path, "trwl_fishery_male_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "female")
ggsave(filename = file.path(fig_path, "trwl_fishery_femlae_mean_LF.png"), width = 8, height = 6)


tag_pred = get_tag_recovery_obs_fitted_values(MLE_report = mle_report, region_key = region_key)
tag_aggregated = tag_pred %>% group_by(release_region, recovery_region) %>% summarise(observed = sum(observed), predicted = sum(predicted))
tag_aggregated$resid = tag_aggregated$observed - tag_aggregated$predicted
tag_aggregated$resid_sign = ifelse(tag_aggregated$resid < 0, "negative", "positive")
tag_aggregated$release_region = factor(tag_aggregated$release_region, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
tag_aggregated$recovery_region = factor(tag_aggregated$recovery_region, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)

gplt = ggplot(tag_aggregated, aes(x = release_region, y = recovery_region, fill = resid)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "steelblue") +
  geom_text(aes(x = release_region, y = recovery_region, label = round(resid,1)), color = "black", size = 4) +
  labs(x = "Recovery", y = "Release", fill = "Aggregated\nresiduals\n (O - E)")
gplt
ggsave(plot = gplt, filename = file.path(fig_path, "tag_aggregated_fit.png"), width = 8, height = 6)


###### more detailed fits
## individual AF's
first_year_set = c(1981, seq(from = 1985, to = 1993, by = 2), 1996:1999)
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "male") +
  ggtitle("Male survey AF") +
  guides(shape = "none", linetype = "none")

plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = first_year_set, sex = "female") +
  ggtitle("Female survey AF") +
  guides(shape = "none", linetype = "none")
second_year_set = 2000:2010
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = second_year_set, sex = "male") +
  ggtitle("Male survey AF") +
  guides(shape = "none", linetype = "none")
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = second_year_set, sex = "female") +
  ggtitle("Female survey AF") +
  guides(shape = "none", linetype = "none")
third_year_set = 2011:2021
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = third_year_set, sex = "male") +
  ggtitle("Male survey AF") +
  guides(shape = "none", linetype = "none") +
  ylim(0,3)
plot_AF(MLE_report = mle_report, region_key = region_key, label = "srv_dom_ll", subset_years = third_year_set, sex = "female") +
  ggtitle("Female survey AF") +
  guides(shape = "none", linetype = "none")

## fishery LFs
# Fixed gear
plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1999, sex = "male") +
  ggtitle("Male fixed LF") +
  guides(shape = "none", linetype = "none")
plot_LF(MLE_report = mle_report, region_key = region_key, label = "fixed", subset_years = 1991:1999, sex = "female") +
  ggtitle("Female fixed LF") +
  guides(shape = "none", linetype = "none")

# Trawl gear
plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 1991:1999, sex = "male") +
  ggtitle("Male trawl LF") +
  guides(shape = "none", linetype = "none") 

plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 1991:1999, sex = "female") +
  ggtitle("Female trawl LF") +
  guides(shape = "none", linetype = "none") 


plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 2000:2010, sex = "male") +
  ggtitle("Male trawl LF") +
  guides(shape = "none", linetype = "none")

plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 2000:2010, sex = "female") +
  ggtitle("Female trawl LF") +
  guides(shape = "none", linetype = "none")

plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 2011:2021, sex = "male") +
  ggtitle("Male trawl LF") +
  guides(shape = "none", linetype = "none")

plot_LF(MLE_report = mle_report, region_key = region_key, label = "trwl", subset_years = 2011:2021, sex = "female") +
  ggtitle("Female trawl LF") +
  guides(shape = "none", linetype = "none") 



#####################
## Find reference points
#####################
proj_data = setup_proj_data(mle_obj = mle_obj, n_proj_years  = 100)
validate_input_data_and_parameters(data = proj_data, parameters = parameters)
proj_pars = mle_spatial$par
proj_obj <- MakeADFun(proj_data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

region_F_ref = find_regional_fref(proj_obj, proj_pars, percent_Bzero = 40, trace = T)
region_F_ref$F_ref
round(region_F_ref$F_ref, 2)

plot_recruitment(MLE_report = region_F_ref$proj_rep, region_key = region_key)
plot_SSB(MLE_report = region_F_ref$proj_rep, region_key = region_key) +
  ylim(0, NA)
plot_SSB(MLE_report = region_F_ref$proj_rep, region_key = region_key, depletion  = T) +
  ylim(0, NA) +
  geom_hline(yintercept = 40, col = "gray60", linetype = "dashed")


## re run with Francis re-weighting
mle_lst = list()
weight_df = NULL
n_reweights = 4
for(i in 1:n_reweights) {
  if(i == 1) {
    ## Re-weight composition observations
    data_weighting = Francis_reweighting(mle_report, region_key)
    # Create factors so we can pivot wider
    data_weighting$length_multipliers$Region = factor(data_weighting$length_multipliers$Region, levels = region_key$area[order(region_key$TMB_ndx)], ordered = T)
    data_weighting$age_multipliers$Region = factor(data_weighting$age_multipliers$Region, levels = region_key$area[order(region_key$TMB_ndx)], ordered = T)
    
    ## pivot_wider for each observations
    fixed_AF_multipliers = data_weighting$age_multipliers %>% filter(label == "fixed") %>% arrange(Region)
    survey_AF_multipliers = data_weighting$age_multipliers %>% filter(label == "srv_dom_ll") %>% arrange(Region)
    fixed_LF_multipliers = data_weighting$length_multipliers %>% filter(label == "fixed") %>% arrange(Region)
    trwl_LF_multipliers = data_weighting$length_multipliers %>% filter(label == "trwl") %>% arrange(Region)
    
    tmp_df = rbind(fixed_AF_multipliers, survey_AF_multipliers, fixed_LF_multipliers, trwl_LF_multipliers)
    tmp_df$iter = i
    weight_df = rbind(weight_df, tmp_df)
    
    # extrapolate and multiple observations by this
    data$obs_fixed_catchatage = sweep(data$obs_fixed_catchatage, MARGIN = c(1,3), FUN = "*", STATS = fixed_AF_multipliers$multiplier)
    data$obs_srv_dom_ll_catchatage = sweep(data$obs_srv_dom_ll_catchatage, MARGIN = c(1,3), FUN = "*", STATS = survey_AF_multipliers$multiplier)
    data$obs_fixed_catchatlgth = sweep(data$obs_fixed_catchatlgth, MARGIN = c(1,3), FUN = "*", STATS = fixed_LF_multipliers$multiplier)
    data$obs_trwl_catchatlgth = sweep(data$obs_trwl_catchatlgth, MARGIN = c(1,3), FUN = "*", STATS = trwl_LF_multipliers$multiplier)
  }
  if(validate_input_data_and_parameters(data, parameters)) {
    obj <- MakeADFun(data, parameters, DLL="SpatialSablefishAssessment_TMBExports")
    mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T)
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(k in 1:2) {
                               g = as.numeric(mle_obj$gr(mle_spatial$par))
                               h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                               mle_spatial$par = mle_spatial$par - solve(h,g)
                               mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                             }
                           , error = function(e){e})
    
    try_improve
    if(inherits(try_improve, "error")) {
      print("non-convergence")
      break;
    }
    mle_report = mle_obj$report(mle_spatial$par)
    mle_lst[[i]] = mle_report
    ## re-weight 
    ## Re-weight composition observations
    data_weighting = Francis_reweighting(mle_report, region_key)
    # Create factors so we can pivot wider
    data_weighting$length_multipliers$Region = factor(data_weighting$length_multipliers$Region, levels = region_key$area[order(region_key$TMB_ndx)], ordered = T)
    data_weighting$age_multipliers$Region = factor(data_weighting$age_multipliers$Region, levels = region_key$area[order(region_key$TMB_ndx)], ordered = T)
    
    ## pivot_wider for each observations
    fixed_AF_multipliers = data_weighting$age_multipliers %>% filter(label == "fixed") %>% arrange(Region)
    survey_AF_multipliers = data_weighting$age_multipliers %>% filter(label == "srv_dom_ll") %>% arrange(Region)
    fixed_LF_multipliers = data_weighting$length_multipliers %>% filter(label == "fixed") %>% arrange(Region)
    trwl_LF_multipliers = data_weighting$length_multipliers %>% filter(label == "trwl") %>% arrange(Region)
    tmp_df = rbind(fixed_AF_multipliers, survey_AF_multipliers, fixed_LF_multipliers, trwl_LF_multipliers)
    tmp_df$iter = i + 1
    weight_df = rbind(weight_df, tmp_df)
    
    # extrapolate and multiple observations by this
    data$obs_fixed_catchatage = sweep(data$obs_fixed_catchatage, MARGIN = c(1,3), FUN = "*", STATS = fixed_AF_multipliers$multiplier)
    data$obs_srv_dom_ll_catchatage = sweep(data$obs_srv_dom_ll_catchatage, MARGIN = c(1,3), FUN = "*", STATS = survey_AF_multipliers$multiplier)
    data$obs_fixed_catchatlgth = sweep(data$obs_fixed_catchatlgth, MARGIN = c(1,3), FUN = "*", STATS = fixed_LF_multipliers$multiplier)
    data$obs_trwl_catchatlgth = sweep(data$obs_trwl_catchatlgth, MARGIN = c(1,3), FUN = "*", STATS = trwl_LF_multipliers$multiplier)
    
  } else {
    print("input error")
    break;
  }
}
## look at the effect of weighting on fitted values
full_index_df = fixed_mean_age_df = survey_mean_age_df = fixed_mean_length_df  = trwl_mean_length_df = NULL
for(i in 1:n_reweights) {
  tmp_ndx_df = get_index(mle_lst[[i]], region_key = region_key)
  tmp_ndx_df$weight = i
  full_index_df = rbind(full_index_df, tmp_ndx_df)
  tmp_srv_mean_age_df = plot_mean_age(MLE_report = mle_lst[[i]], region_key = region_key, label = "srv_dom_ll", sex = "both")$data
  tmp_srv_mean_age_df$weight = i
  survey_mean_age_df = rbind(survey_mean_age_df, tmp_srv_mean_age_df)
  tmp_fixed_mean_age_df = plot_mean_age(MLE_report = mle_lst[[i]], region_key = region_key, label = "fixed", sex = "both")$data
  tmp_fixed_mean_age_df$weight = i
  fixed_mean_age_df = rbind(fixed_mean_age_df, tmp_fixed_mean_age_df)
  tmp_fixed_mean_len_df = plot_mean_length(MLE_report = mle_lst[[i]], region_key = region_key, label = "fixed", sex = "both")$data
  tmp_fixed_mean_len_df$weight = i
  fixed_mean_length_df = rbind(fixed_mean_length_df, tmp_fixed_mean_len_df)
  tmp_trwl_mean_len_df = plot_mean_length(MLE_report = mle_lst[[i]], region_key = region_key, label = "trwl", sex = "both")$data
  tmp_trwl_mean_len_df$weight = i
  trwl_mean_length_df = rbind(trwl_mean_length_df, tmp_trwl_mean_len_df)
}

## re-estimate


## simdata
simdata = mle_obj$simulate(par = mle_spatial$par, complete = T)
OM_pars = parameters
save(simdata, OM_pars, file = file.path(fig_path, "simulated_data.RData"))
## 
years[43]
profile_rec_2019 <- profile_param(parameters = parameters, mle_obj = mle_obj, na_map = map_fixed_pars, profile_param_label = "ln_rec_dev", 
                                  element = matrix(c(1,43), nrow = 1), profile_values = c(-2, -1,0,2,4,6,8),  same_element = NULL)

profile_param_label = "ln_rec_dev"
profile_values = c(-2, -1,0,2,4,6,8)
element = matrix(c(1,43), nrow = 1)
na_map = map_fixed_pars
same_element = NULL
