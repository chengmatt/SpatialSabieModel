#'
#' Single area model using my data methods 
#' Differs from current assessemnt by starting in 1977 instead of 1960
#' We also don't include the early observations and
#' we have restructured the composition observations to be sexually disaggregared.

source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
# options(warn=0) ## if warnings are being set to errors
scenario = "1977_start"
AF_direct_ageing = T
max_N_eff = 500 ## during the initial input calculations some of the input sample sizes were huge, which I don't want effecting the model
min_N_eff = 50

N_eff_multiplier = 0.4
fig_path = file.path(DIR$app_1A, scenario)
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}

## read in observational datasets
fixed_gear_AF_alk_pooled = readRDS(file = file.path(DIR$data1A, "Observer_ALK_AF_w_eff.RDS"))
fixed_gear_AF_direct = readRDS(file = file.path(DIR$data1A, "Observer_direct_AF_w_eff.RDS"))
trawl_observer_LF = readRDS(file = file.path(DIR$data1A, "Observer_trawl_LF_w_eff.RDS"))
fixed_observer_LF = readRDS(file = file.path(DIR$data1A, "Observer_fixed_LF_w_eff.RDS"))
#survery_AF_direct = readRDS(file = file.path(DIR$data1A, "Survey_direct_AF_w_eff.RDS"))
#survery_AF_alk_pooled = readRDS(file = file.path(DIR$data1A, "Survey_ALK_AF_w_eff.RDS"))
survery_AF_direct = readRDS(file = file.path(DIR$data1A, "Survey_direct_ByCountry_AF_w_eff.RDS"))
survery_AF_alk_pooled = readRDS(file = file.path(DIR$data1A, "Survey_ALK_ByCountry_AF_w_eff.RDS"))


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
                                                       case_when(
                                                         area_lab == "AI" ~ "SingleArea",
                                                         area_lab == "BS" ~ "SingleArea",
                                                         area_lab == "WGOA" ~ "SingleArea",
                                                         area_lab == "EGOA" ~ "SingleArea",
                                                         area_lab == "CGOA" ~ "SingleArea"))
## sum estimates over all years
design_survey_index = design_survey_index %>% group_by(Year, area_lab) %>% summarise(sum_estimates = sum(sum_estimates), sum_var = sum(sum_var), se = sqrt(sum_var), LCI = sum_estimates - 2*se, UCI = sum_estimates + 2*se)

ggplot(design_survey_index, aes(Year, sum_estimates, col = area_lab, linetype = area_lab, fill = area_lab, alpha = 0.3)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI)) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Relative index (RPN)", linetype = "", col = "", fill = "") +
  guides(alpha = "none") +
  facet_wrap(~area_lab, scales = "free_y") +
  theme_bw()
ggsave(filename = file.path(DIR$input_figs_1A, "indicies.png"), width = 7, height = 7)

srv_jap_df = data.frame(years = sab_inputs$srv_jap_ll_bio_yrs, obs =  sab_inputs$srv_jap_ll_bio_obs, se = sab_inputs$srv_jap_ll_bio_se, label = "Japanese LL Survey")
srv_dom_df = data.frame(years = sab_inputs$srv_dom_ll_bio_yrs, obs =  sab_inputs$srv_dom_ll_obs, se = sab_inputs$srv_dom_ll_se, label = "Domestic LL Survey")
srv_nmfs_trwl_df = data.frame(years = sab_inputs$srv_nmfs_trwl_bio_yrs, obs =  sab_inputs$srv_nmfs_trwl_obs, se = sab_inputs$srv_nmfs_trwl_se, label = "NMFS Trawl Survey")
japanese_fishery_cpue_df = data.frame(years = sab_inputs$srv_jap_fishery_ll_bio_yrs, obs =  sab_inputs$srv_jap_fishery_ll_obs, se = sab_inputs$srv_jap_fishery_ll_se, label = "Japanese fishery CPUE")
full_ass_ndx = rbind(srv_jap_df,srv_dom_df, srv_nmfs_trwl_df, japanese_fishery_cpue_df)
ggplot(full_ass_ndx, aes(x = years, y = obs, col = label, linetype = label)) +
  geom_point(size = 1) +
  geom_line(linewidth = 1.1) +
  theme_bw()
## Tag data
tag_recovery_df = readRDS(file = file.path(DIR$data1A, "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = file.path(DIR$data1A, "Tag_release_summarised.RDS"))

## read in Catch
full_catch_df = readRDS(file = file.path(DIR$data1A, "Catch_by_year_area_gear.RDS"))

## bring ADMB inputs to help configure initial model
sab_curr <- dget(file.path(DIR$admb, "tem.rdat")) 
sab_ctl <- digest_sab_ctl(sablefish_control_filename = file.path(DIR$admb, "tem.ctl")) 
sab_inputs <- digest_sab_input(sablefish_input_filename = file.path(DIR$admb, "tem_2022_na_wh.dat")) 
sab_par <- readLines(file.path(DIR$admb, "tem.par")) 
sab_rep <-readLines(file.path(DIR$admb, "sable.rep"))
tem_rep <-readLines(file.path(DIR$admb, "tem.rep"))
admb_out <-readLines(file.path(DIR$admb, "SABLE_SARA.dat"))
## Model dimensions
admb_years = sab_inputs$yrs[1]:sab_inputs$yrs[2]
years = 1977:2021
min_age = 2
max_age = 31
ages = min_age:max_age
region_key = data.frame(area = c("SingleArea"), TMB_ndx = c(0))
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
survey_labels = c("Early Japanese Fishery", "Japanese LL Survey","Domestic LL survey" )#, "NMFS Trawl Survey")
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
data$n_init_rec_devs =  sab_inputs$yrs[1] - sab_ctl$rec_styr - 2## number of non-equilibrium age-structure
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
data$male_age_length_transition[,,as.character(min(years):1994)] = sab_inputs$male_length_at_age_60_96
data$female_age_length_transition[,,as.character(min(years):1994)] = sab_inputs$female_length_at_age_60_96
data$male_age_length_transition[,,as.character(1995:max(years))] = sab_inputs$male_length_at_age_97_22
data$female_age_length_transition[,,as.character(1995:max(years))] = sab_inputs$female_length_at_age_97_22


## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
data$movement_matrix = data$fixed_movement_matrix = matrix(0, nrow = n_regions, ncol = n_regions);
diag(data$movement_matrix) = 0.9
diag(data$fixed_movement_matrix) = 1
data$movement_matrix = data$movement_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1)
data$movement_matrix = sweep(data$movement_matrix, 1, STATS = rowSums(data$movement_matrix), "/")


# renormalise
data$spawning_time_proportion = rep(0, n_projyears)
data$sigma_R = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("sigr:",sab_par)+1],split=" "))))
data$SrType = 3
data$apply_fixed_movement = 0; ##
#'
#' Fishing inputs
#' 
data$F_method = 0
data$F_max = 3
data$F_iterations = 4
## 
data$prop_F_hist = 0.1
# drop years outside model years
full_catch_df = full_catch_df %>% filter(year %in% years)
full_catch_df$year_f = factor(full_catch_df$year, levels = year_lvls, ordered = T)
full_catch_df$region_f = factor(full_catch_df$area, levels = reg_lvls, ordered = T)
## convert catch to metric tonnes
full_catch_df$catch_mt = full_catch_df$catch / 1000
trwl_catch = full_catch_df %>% filter(fmp_gear == "TRW") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)
fixed_catch = full_catch_df %>% filter(fmp_gear == "HAL") %>% ungroup() %>% dplyr::select(catch_mt, year_f, region_f) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = catch_mt)

## replace 0 values with low values
trwl_catch = data.frame(sab_inputs$trwl_Catch)
fixed_catch = data.frame(sab_inputs$ll_Catch)

fixed_catch = fixed_catch %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))
trwl_catch = trwl_catch %>% mutate(across(.fns = ~replace(., . ==  0 , 0.001)))


data$trwl_fishery_catch = matrix(trwl_catch$sab_inputs.trwl_Catch[sab_inputs$yrs[1]:sab_inputs$yrs[2] %in% years], nrow = n_regions, ncol = n_projyears)
data$fixed_fishery_catch = matrix(fixed_catch$sab_inputs.ll_Catch[sab_inputs$yrs[1]:sab_inputs$yrs[2] %in% years], nrow = n_regions, ncol = n_projyears)

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
data$srv_sel_type = matrix(0, nrow = 1, ncol = data$n_surveys)
#data$srv_sel_type[4] = 4 ## Trawl survey has exponential decay
data$srv_sel_by_year_indicator = matrix(0, nrow = n_projyears, ncol = data$n_surveys)
#data$srv_sel_by_year_indicator[projyears %in% 1988:1994, 1] = 1
#data$srv_sel_by_year_indicator[projyears > 1994, 1] = 2


#'
#'
#' Tag release stuff
#'
tag_release_years = max(years)#, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016)
data$tag_release_event_this_year = rep(0, n_years) ## no tag releases
data$male_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))
data$female_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))

#paste(data$years[which(data$tag_release_event_this_year == 1)], collapse = ", ")

data$n_years_to_retain_tagged_cohorts_for = 1
data$initial_tag_induced_mortality = rep(0, length(tag_release_years))
data$annual_tag_shedding_rate = 0.00

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

data$srv_catchatage_indicator[1,data$years %in% US_survey_years,3] = 1
data$srv_catchatage_indicator[1,data$years %in% jap_survey_years,2] = 1

data$obs_srv_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years, data$n_surveys))
obs_years = unique(survery_AF$year)
obs_reg = unique(survery_AF$region)
## populate container
countries = unique(survery_AF$country)
for(c_ndx in 1:length(countries)) { ## go backwards becuase its starts with Japan
  for(y_ndx in 1:length(obs_years)) {
    this_year_ndx = which(years %in% obs_years[y_ndx])
    for(r_ndx in 1:length(obs_reg)) {
      this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
      this_df = (survery_AF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx], country == countries[c_ndx])) %>% arrange(sex_age)
      if(nrow(this_df) > 0)
        data$obs_srv_catchatage[,this_region_ndx, this_year_ndx, c_ndx + 1] = this_df$P * this_df$eff_N
    }
  }
}
data$srv_catchatage_covar_structure = 0
data$srv_catchatage_comp_likelihood = rep(0, data$n_surveys)

### Survey LL index
data$srv_bio_indicator = array(0, dim = c(n_regions, n_years,data$n_surveys))
data$srv_bio_indicator[1,data$years %in% srv_dom_df$years, 3] = 1
data$srv_bio_indicator[1,data$years %in% srv_jap_df$years, 2] = 1
data$srv_bio_indicator[1,data$years %in% japanese_fishery_cpue_df$years, 1] = 1

data$obs_srv_bio = array(0, dim = c(n_regions, n_years, data$n_surveys))
data$obs_srv_bio[1,data$years %in% srv_dom_df$years, 3] = srv_dom_df$obs[srv_dom_df$years %in% data$years]
data$obs_srv_bio[1,data$years %in% srv_jap_df$years, 2] = srv_jap_df$obs[srv_jap_df$years %in% data$years]
data$obs_srv_bio[1,data$years %in% japanese_fishery_cpue_df$years, 1] = japanese_fishery_cpue_df$obs[japanese_fishery_cpue_df$years %in% data$years]

data$obs_srv_se = array(0.0, dim = c(n_regions, n_years,data$n_surveys))
data$obs_srv_se[1,data$years %in% srv_dom_df$years, 3] = srv_dom_df$se[srv_dom_df$years %in% data$years]
data$obs_srv_se[1,data$years %in% srv_jap_df$years, 2] = srv_jap_df$se[srv_jap_df$years %in% data$years]
data$obs_srv_se[1,data$years %in% japanese_fishery_cpue_df$years, 1] = japanese_fishery_cpue_df$se[japanese_fishery_cpue_df$years %in% data$years]

data$srv_bio_likelihood = rep(0, data$n_surveys)
data$srv_obs_is_abundance = rep(0, data$n_surveys)
data$srv_q_by_year_indicator = matrix(0, nrow = n_years, ncol = data$n_surveys)
#data$srv_q_by_year_indicator[projyears %in% 1988:1994, 1] = 1
#data$srv_q_by_year_indicator[projyears > 1994, 1] = 2
data$srv_q_transformation = rep(0, data$n_surveys)
data$q_is_nuisance = rep(0, data$n_surveys)

###
# Tag-recoveries
# This can take a bit of time to get your head around!
###
include_tag_recoveries = T
include_zero_tag_recovery_events = T
tag_recovery_years = 1978:2020

data$tag_recovery_indicator_by_year = rep(0, n_years) ## no tag releases
data$obs_tag_recovery = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))
data$tag_recovery_indicator = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))

## are there observatiosn for the plus group
if(sum(data$tag_recovery_indicator[1:data$n_regions,,]) != 0)
  cat("found recoveries in the first year of release\n")
if(sum(data$tag_recovery_indicator[(dim(data$tag_recovery_indicator)[1] - data$n_regions + 1):dim(data$tag_recovery_indicator)[1],,]) != 0)
  cat("found recoveries in the plus group\n")

cat("The number of tag recovery observations ", sum(data$tag_recovery_indicator), "\n")

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
parameters$ln_mean_rec = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_mean_rec:",sab_par)+1],split=" "))))
## Fishery selectivity
parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_trwl_sel_pars = array(0, dim = c(length(unique(data$trwl_sel_by_year_indicator )), 2, 2))
parameters$ln_srv_sel_pars = array(0, dim = c(max(apply(data$srv_sel_by_year_indicator, 2,FUN = function(x){length(unique(x))})), 2, 2, data$n_surveys))

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
if(dim(parameters$ln_srv_sel_pars)[1] == 1) {
  parameters$ln_srv_sel_pars[1,1,1,] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1,] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2,] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2,] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
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
} else if(dim(parameters$ln_srv_sel_pars)[1] == 4) {
  parameters$ln_srv_sel_pars[1,1,1] =parameters$ln_srv_sel_pars[2,1,1] = parameters$ln_srv_sel_pars[3,1,1] = parameters$ln_srv_sel_pars[4,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,1] = parameters$ln_srv_sel_pars[2,2,1] = parameters$ln_srv_sel_pars[3,2,1] = parameters$ln_srv_sel_pars[4,2,1] =suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,1,2] = parameters$ln_srv_sel_pars[2,1,2] = parameters$ln_srv_sel_pars[3,1,2] = parameters$ln_srv_sel_pars[4,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
  parameters$ln_srv_sel_pars[1,2,2] = parameters$ln_srv_sel_pars[2,2,2] = parameters$ln_srv_sel_pars[3,2,2] = parameters$ln_srv_sel_pars[4,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))
}

## movement pars
if(n_regions > 1) {
  parameters$transformed_movement_pars = matrix(NA, nrow = n_regions - 1, ncol = n_regions)
  for(i in 1:n_regions)
    parameters$transformed_movement_pars[,i] = simplex(data$movement_matrix[i,])
} else {
  parameters$transformed_movement_pars = matrix(1, nrow = 1, ncol = 1)
}
parameters$ln_fixed_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish1:",sab_par)+1],split=" "))))
admb_fixed_devs = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish1:",sab_par)+1],split=" "))))[-1]
parameters$ln_fixed_F_devs = matrix(admb_fixed_devs[admb_years %in% years], nrow = n_regions)

parameters$ln_trwl_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish3:",sab_par)+1],split=" "))))
admb_trwl_devs = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish3:",sab_par)+1],split=" "))))[-1]
parameters$ln_trwl_F_devs = matrix(admb_trwl_devs[admb_years %in% years], nrow = n_regions)

parameters$ln_init_F_avg = parameters$ln_fixed_F_avg
parameters$trans_srv_q = array(log(1000), dim =c(n_regions,max(apply(data$srv_q_by_year_indicator, 2,FUN = function(x){length(unique(x))})), data$n_surveys))

parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years)


#exp(ln_mean_rec + ln_rec_dev(year_ndx) + sigma_R_sq/2
#parameters$ln_rec_dev = rep(0 - 0.5 * data$sigma_R * data$sigma_R, n_projyears)
parameters$ln_catch_sd = log(0.02)
if(ncol(parameters$trans_rec_dev) != n_years)
  stop("input error with trans_rec_dev")

if(data$n_init_rec_devs == 0) {
  parameters$ln_init_rec_dev = 0#rev(admb_ln_rec_devs[1:data$n_init_rec_devs]) ## note the reverse fucntion here. see source code for why we have to do this
} else {
  parameters$ln_init_rec_dev = rep(0, data$n_init_rec_devs)
}
parameters$logistic_tag_reporting_rate = matrix(logit(0.50), nrow = data$n_regions, ncol = max(sum(data$tag_recovery_indicator_by_year),1))

parameters$ln_tag_phi = log(0.5)
parameters$ln_sigma_R = log(data$sigma_R)
parameters$ln_sigma_init_devs = log(1.2)

parameters$trans_trwl_catchatlgth_error = log(1)
parameters$trans_fixed_catchatlgth_error = log(1)
parameters$trans_fixed_catchatage_error = log(1)
parameters$trans_srv_catchatage_error = rep(log(1), data$n_surveys)
parameters$trans_SR_pars = rep(log(0.8),1)
parameters$logistic_prop_recruit_male = rep(logit(0.5), length(data$years))

######
## Check inputs
######
data$model = "TagIntegrated"
data$F_method = 1 ## estimate free F's
data$apply_fixed_movement = 1
data$do_recruits_move = 0

data$prop_F_hist = 1.0
data$tag_likelihood = 1
data$evaluate_tag_likelihood = 1

#data$n_init_rec_devs = 0

validate_input_data_and_parameters(data, parameters)

## save objects to speed up time in debugging
#
if(FALSE) {
  data = readRDS(file.path(fig_path, "data.RDS"))
  parameters = readRDS(file.path(fig_path, "parameters.RDS"))
}

## plot some input info
plot_input_observations(data, region_key = region_key) + theme(axis.text = element_text(size = 14))
ggsave(filename = file.path(fig_path, "Observation_Frequency.png"), width = 7, height = 7)

plot_input_timeblocks(data)
ggsave(filename = file.path(fig_path, "time_blocks.png"), width = 7, height = 6)


plot_input_catches(data, region_key = region_key) + theme(axis.text = element_text(size = 14))
ggsave(filename = file.path(fig_path, "InputCatches.png"), width = 7, height = 7)

#####################
## try and estimate 
#####################
## some pars to fix
#na_map = fix_pars(par_list = parameters, pars_to_exclude = "trans_srv_dom_ll_q")

map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = F,
                                   srv_sel_second_param_shared_by_sex = T,
                                   fixed_sel_first_shared_by_sex  = F,
                                   fixed_sel_second_shared_by_sex   = T,
                                   trwl_sel_first_shared_by_sex  = F,
                                   trwl_sel_second_shared_by_sex  = F,
                                   recruit_dev_years_not_to_estimate = NULL, 
                                   srv_q_spatial = F,
                                   tag_reporting_rate = "off",
                                   est_init_F = T,
                                   est_catch_sd = F,
                                   est_movement = F,
                                   est_prop_male_recruit = "off"
                                   #est_prop_male_recruit = seq(from = 1990, to = 2020, by = 10)
)

## Make AD Fun
mle_obj <- MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T, silent = T)

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
sim_obs = simulate_observations(obj = mle_obj, n_sims = 200, sd_report = sd_report, include_param_uncertainty = F, region_key = region_key)
## save this so we don't have to keep running it
saveRDS(sim_obs, file.path(fig_path, "sim_obs.RDS"))

plot_comp_sample_size(mle_report, data, region_key, F)
# profile R0
ln_R0_profile_pars = seq(from = 13, to = 16, by = 0.1)
profile_ln_R0 = profile_param(mle_obj = mle_obj, parameters = mle_param_list, na_map = map_fixed_pars, profile_param_label = "ln_mean_rec", profile_values = ln_R0_profile_pars)
saveRDS(profile_ln_R0, file.path(fig_path, "profile_ln_R0.RDS"))


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
ggsave(filename = file.path(fig_path, "Init_numbers_at_age.png"), width = 5, height = 5)

# plot selectivities
plot_selectivities(mle_report)
ggsave(filename = file.path(fig_path, "Est_selectivity.png"), width = 7, height = 6)

# starting values
plot_selectivities(start_report)

## SSBs
ssb_plt = plot_SSB(mle_report, region_key = region_key)
ggsave(filename = file.path(fig_path, "regional_ssb.png"), width = 6, height = 6)
plot_SSB(mle_report, region_key = region_key, depletion = T)
ggsave(filename = file.path(fig_path, "depletion_ssb.png"), width = 6, height = 6)


## sum SSB over all regions
ggplot(ssb_plt$data %>% group_by(Year) %>% summarise(total_ssb = sum(SSB)), aes(x = Year, y = total_ssb)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA)
ggsave(filename = file.path(fig_path, "global_ssb.png"), width = 6, height = 6)

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
ggsave(filename = file.path(fig_path, "survey_male_mean_AF.png"), width = 6, height = 5)
plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "srv_dom_ll",sex = "female")
ggsave(filename = file.path(fig_path, "survey_female_mean_AF.png"), width = 6, height = 5)

plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male")
ggsave(filename = file.path(fig_path, "fixed_fishery_male_mean_AF.png"), width = 6, height =5)
plot_mean_age(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female")
ggsave(filename = file.path(fig_path, "fixed_fishery_femlae_mean_AF.png"), width = 6, height = 5)

plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "male")
ggsave(filename = file.path(fig_path, "fixed_fishery_male_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "fixed",sex = "female")
ggsave(filename = file.path(fig_path, "fixed_fishery_femlae_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "male")
ggsave(filename = file.path(fig_path, "trwl_fishery_male_mean_LF.png"), width = 8, height = 6)
plot_mean_length(MLE_report = mle_report, region_key = region_key,label = "trwl",sex = "female")
ggsave(filename = file.path(fig_path, "trwl_fishery_femlae_mean_LF.png"), width = 8, height = 6)


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
  guides(shape = "none", linetype = "none") 
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

########################
## Compare TMB output with ADMB
########################
## Extact ADMB output
mean_recruit<-as.numeric(unlist(strsplit(sab_rep[grep("Mean_Recruitment",sab_rep)+1],split=" "))) #Get mean recruitment used for B40 calc from report file
SSB <- as.numeric(unlist(strsplit(sab_rep[grep("SpBiom",sab_rep)],split=" ")))[-c(1,2)] #Get mean recruitment used for B40 calc from report file
catch_fish1 <- as.numeric(unlist(strsplit(admb_out[grep("#CATCHFISH1 - catch by year in kilo tons",admb_out) +1],split=" ")))[-1] #Get mean recruitment used for B40 calc from report file
catch_fish3 <- as.numeric(unlist(strsplit(admb_out[grep("#CATCHFISH3 - catch by year in kilo tons",admb_out) +1],split=" ")))[-1] #Get mean recruitment used for B40 calc from report file
N_age_f<-as.numeric(unlist(strsplit(admb_out[grep("N_AT_AGE",admb_out)+1],split=" "))) #Get numbers at age
N_age_m<-as.numeric(unlist(strsplit(admb_out[grep("N_AT_AGE",admb_out)+2],split=" "))) #Get numbers at age
N_length_f<-as.numeric(unlist(strsplit(sab_rep[(grep("Numbers at Length Females",sab_rep)+1):(grep("Numbers at Length Females",sab_rep)+n_years)],split=" "))) #Get numbers at age
N_length_m<-as.numeric(unlist(strsplit(sab_rep[(grep("Numbers at Length  Males",sab_rep)+1):(grep("Numbers at Length  Males",sab_rep)+n_years)],split=" "))) #Get numbers at age
catchatage_ll_f<-as.numeric(unlist(strsplit(admb_out[(grep("#CatchAtAge_FISH1_f",admb_out)+1):(grep("#CatchAtAge_FISH1_f",admb_out)+n_years)],split=" "))) #Get numbers at age
catchatage_ll_m<-as.numeric(unlist(strsplit(admb_out[(grep("#CatchAtAge_FISH1_m",admb_out)+1):(grep("#CatchAtAge_FISH1_m",admb_out)+n_years)],split=" "))) #Get numbers at age
catchatage_trwl_f<-as.numeric(unlist(strsplit(admb_out[(grep("#CatchAtAge_FISH3_f",admb_out)+1):(grep("#CatchAtAge_FISH3_f",admb_out)+n_years)],split=" "))) #Get numbers at age
catchatage_trwl_m<-as.numeric(unlist(strsplit(admb_out[(grep("#CatchAtAge_FISH3_m",admb_out)+1):(grep("#CatchAtAge_FISH3_m",admb_out)+n_years)],split=" "))) #Get numbers at age

SSB_rec <- (unlist(strsplit(sab_rep[grep("Year SSB SRR Recr",sab_rep):(grep("Year SSB SRR Recr",sab_rep) + length(admb_years))],split=" "))) 
SSB_rec_mat = matrix(SSB_rec, byrow = T, ncol = 4)
header = SSB_rec_mat[1,]
SSB_rec_mat = SSB_rec_mat[-1,]
colnames(SSB_rec_mat) = header
class(SSB_rec_mat) = "numeric"

Fatage_ll_f<-as.numeric(unlist(strsplit(admb_out[(grep("#FAtAge_FISH1_f",admb_out)+1):(grep("#FAtAge_FISH1_f",admb_out)+n_years)],split=" "))) #Get numbers at age
Fatage_ll_m<-as.numeric(unlist(strsplit(admb_out[(grep("#FAtAge_FISH1_m",admb_out)+1):(grep("#FAtAge_FISH1_m",admb_out)+n_years)],split=" "))) #Get numbers at age
Fatage_trwl_f<-as.numeric(unlist(strsplit(admb_out[(grep("#FAtAge_FISH3_f",admb_out)+1):(grep("#FAtAge_FISH3_f",admb_out)+n_years)],split=" "))) #Get numbers at age
Fatage_trwl_m<-as.numeric(unlist(strsplit(admb_out[(grep("#FAtAge_FISH3_m",admb_out)+1):(grep("#FAtAge_FISH3_m",admb_out)+n_years)],split=" "))) #Get numbers at age
F_ll_yr <- as.numeric(unlist(strsplit(admb_out[(grep("#Fmort_fish1",admb_out)+1)],split=" ")))[-1] #Get numbers at age
F_trwl_yr <- as.numeric(unlist(strsplit(admb_out[(grep("#Fmort_fish3",admb_out)+1)],split=" ")))[-1] #Get numbers at age

sel_ll_f<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+1)],split=" "))) #Get numbers at age
sel_ll_m<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+2)],split=" "))) #Get numbers at age
sel4_ll_f<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+3)],split=" "))) #Get numbers at age
sel4_ll_m<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+4)],split=" "))) #Get numbers at age
sel5_ll_f<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+5)],split=" "))) #Get numbers at age
sel5_ll_m<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated fixed gear fishery selectivity",admb_out)+6)],split=" "))) #Get numbers at age
sel_trwl_f<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated trawl gear fishery selectivity",admb_out)+1)],split=" "))) #Get numbers at age
sel_trwl_m<-as.numeric(unlist(strsplit(admb_out[(grep("Estimated trawl gear fishery selectivity",admb_out)+2)],split=" "))) #Get numbers at age

sel_srv_dom_ll_1_f<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey domestic_LL gear fishery selectivity",admb_out)+1)],split=" "))) #Get numbers at age
sel_srv_dom_ll_1_m<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey domestic_LL gear fishery selectivity",admb_out)+2)],split=" "))) #Get numbers at age
sel_srv_dom_ll_10_f<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey domestic_LL gear fishery selectivity",admb_out)+3)],split=" "))) #Get numbers at age
sel_srv_dom_ll_10_m<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey domestic_LL gear fishery selectivity",admb_out)+4)],split=" "))) #Get numbers at age

sel_srv_jap_ll_2_f<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey japanese_LL gear fishery selectivity",admb_out)+1)],split=" "))) #Get numbers at age
sel_srv_jap_ll_2_m<-as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey japanese_LL gear fishery selectivity",admb_out)+2)],split=" "))) #Get numbers at age

sel_srv_nmfs_trwl_f = as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey Trawl GOA gear fishery selectivit",admb_out)+1)],split=" "))) #Get numbers at age
sel_srv_nmfs_trwl_m = as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated survey Trawl GOA gear fishery selectivit",admb_out)+2)],split=" "))) #Get numbers at age

sel_srv_jap_fishery_ll_m = as.numeric(unlist(strsplit(admb_out[(grep("#SELECTIVITY - Estimated historic japanese_LL  fishery selectivity",admb_out)+1)],split=" "))) #Get numbers at age

## initial numbers at age
init_N_age_f<-as.numeric(unlist(strsplit(admb_out[grep("init_Numbers_AT_AGE",admb_out)+1],split=" "))) #Get numbers at age
init_N_age_m<-as.numeric(unlist(strsplit(admb_out[grep("init_Numbers_AT_AGE",admb_out)+2],split=" "))) #Get numbers at age
mdelta<-as.numeric(unlist(strsplit(admb_out[grep("mdelta",admb_out)+1],split=" "))) #Get numbers at age

# drop NA's 
init_N_age_f = subset(init_N_age_f, subset = !is.na(init_N_age_f))
init_N_age_m = subset(init_N_age_m, subset = !is.na(init_N_age_m))
N_age_f = subset(N_age_f, subset = !is.na(N_age_f))
N_age_m = subset(N_age_m, subset = !is.na(N_age_m))
catchatage_ll_m = subset(catchatage_ll_m, subset = !is.na(catchatage_ll_m))
catchatage_ll_f = subset(catchatage_ll_f, subset = !is.na(catchatage_ll_f))
catchatage_trwl_m = subset(catchatage_trwl_m, subset = !is.na(catchatage_trwl_m))
catchatage_trwl_f = subset(catchatage_trwl_f, subset = !is.na(catchatage_trwl_f))


Fatage_ll_m = subset(Fatage_ll_m, subset = !is.na(Fatage_ll_m))
Fatage_ll_f = subset(Fatage_ll_f, subset = !is.na(Fatage_ll_f))
Fatage_trwl_m = subset(Fatage_trwl_m, subset = !is.na(Fatage_trwl_m))
Fatage_trwl_f = subset(Fatage_trwl_f, subset = !is.na(Fatage_trwl_f))

N_length_f = subset(N_length_f, subset = !is.na(N_length_f))
N_length_m = subset(N_length_m, subset = !is.na(N_length_m))
sel_ll_f = subset(sel_ll_f, subset = !is.na(sel_ll_f))
sel_ll_m = subset(sel_ll_m, subset = !is.na(sel_ll_m))
sel4_ll_f = subset(sel4_ll_f, subset = !is.na(sel4_ll_f))
sel4_ll_m = subset(sel4_ll_m, subset = !is.na(sel4_ll_m))
sel5_ll_f = subset(sel5_ll_f, subset = !is.na(sel5_ll_f))
sel5_ll_m = subset(sel5_ll_m, subset = !is.na(sel5_ll_m))
sel_trwl_f = subset(sel_trwl_f, subset = !is.na(sel_trwl_f))
sel_trwl_m = subset(sel_trwl_m, subset = !is.na(sel_trwl_m))

sel_srv_dom_ll_1_f = subset(sel_srv_dom_ll_1_f, subset = !is.na(sel_srv_dom_ll_1_f))
sel_srv_dom_ll_1_m = subset(sel_srv_dom_ll_1_m, subset = !is.na(sel_srv_dom_ll_1_m))
sel_srv_dom_ll_10_f = subset(sel_srv_dom_ll_10_f, subset = !is.na(sel_srv_dom_ll_10_f))
sel_srv_dom_ll_10_m = subset(sel_srv_dom_ll_10_m, subset = !is.na(sel_srv_dom_ll_10_m))

sel_srv_jap_ll_2_f = subset(sel_srv_jap_ll_2_f, subset = !is.na(sel_srv_jap_ll_2_f))
sel_srv_jap_ll_2_m = subset(sel_srv_jap_ll_2_m, subset = !is.na(sel_srv_jap_ll_2_m))

sel_srv_nmfs_trwl_f = subset(sel_srv_nmfs_trwl_f, subset = !is.na(sel_srv_nmfs_trwl_f))
sel_srv_nmfs_trwl_m = subset(sel_srv_nmfs_trwl_m, subset = !is.na(sel_srv_nmfs_trwl_m))
sel_srv_jap_fishery_ll = subset(sel_srv_jap_fishery_ll_m, subset = !is.na(sel_srv_jap_fishery_ll_m))
# turn to matrix
N_age_f_mat = matrix(N_age_f, byrow = T, nrow = n_years, ncol = n_ages)
N_age_m_mat = matrix(N_age_m, byrow = T, nrow = n_years, ncol = n_ages)

catchatage_ll_m_mat = matrix(catchatage_ll_m, byrow = T, nrow = n_years, ncol = n_ages)
catchatage_ll_f_mat = matrix(catchatage_ll_f, byrow = T, nrow = n_years, ncol = n_ages)
catchatage_trwl_m_mat = matrix(catchatage_trwl_m, byrow = T, nrow = n_years, ncol = n_ages)
catchatage_trwl_f_mat = matrix(catchatage_trwl_f, byrow = T, nrow = n_years, ncol = n_ages)

Fatage_ll_m_mat = matrix(Fatage_ll_m, byrow = T, nrow = n_years, ncol = n_ages)
Fatage_ll_f_mat = matrix(Fatage_ll_f, byrow = T, nrow = n_years, ncol = n_ages)
Fatage_trwl_m_mat = matrix(Fatage_trwl_m, byrow = T, nrow = n_years, ncol = n_ages)
Fatage_trwl_f_mat = matrix(Fatage_trwl_f, byrow = T, nrow = n_years, ncol = n_ages)

N_length_f_mat = matrix(N_length_f, byrow = T, nrow = n_years, ncol = n_length_bins + 1)
N_length_m_mat = matrix(N_length_m, byrow = T, nrow = n_years, ncol = n_length_bins + 1)
## drop leading column (its the year)
N_length_f_mat = N_length_f_mat[,-1]
N_length_m_mat = N_length_m_mat[,-1]

## plot SSB
png(filename = file.path(fig_path, "compare_SSB.png"), width = 6, height = 6, units ="in", res = 250)
plot(admb_years, SSB, lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "SSB", ylim = c(0,400))
lines(data$years, mle_report$SSB, lwd = 3, col = "black", lty = 2)
legend('topright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd = 2)
dev.off()

## plot Annual F's
png(filename = file.path(fig_path, "compare_Fs.png"), width = 8, height = 8, units ="in", res = 250)
par(mfrow = c(2,1))
plot(admb_years, F_ll_yr, lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "F", main = "Longline", ylim = c(0,0.15))
lines(data$years, mle_report$annual_F_fixed, lwd = 3, col = "black", lty = 2)
plot(admb_years, F_trwl_yr, lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "F", main = "Trawl")
lines(data$years, mle_report$annual_F_trwl, lwd = 3, col = "black", lty = 2)
legend('topright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

## plot Annual recruitment
png(filename = file.path(fig_path, "compare_Recruitment.png"), width = 6, height = 6, units ="in", res = 250)
par(mfrow = c(1,1))
plot(admb_years, SSB_rec_mat[,"Recr"], lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "Recruits", ylim = c(0,100))
lines(data$years, mle_report$recruitment_yr[,1], lwd = 3, col = "black", lty = 2)
legend('topleft', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

## Compare selectivties
# first LL selectivity
png(filename = file.path(fig_path, "compare_sel_LL_1_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel_ll_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_fixed_f[,1], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel_ll_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Male")
lines(data$ages, mle_report$sel_fixed_m[,1], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

png(filename = file.path(fig_path, "sel_LL_4_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel4_ll_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_fixed_f[,2], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel4_ll_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Male")
lines(data$ages, mle_report$sel_fixed_m[,2], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

png(filename = file.path(fig_path, "sel_LL_5_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel5_ll_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_fixed_f[,3], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel5_ll_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Male")
lines(data$ages, mle_report$sel_fixed_m[,3], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

png(filename = file.path(fig_path, "sel_trwl_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel_trwl_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_trwl_f[,1], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel_trwl_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Mmale")
lines(data$ages, mle_report$sel_trwl_m[,1], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

# Survey dom LL selectivity
png(filename = file.path(fig_path, "sel_srv_dom_LL_1_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel_srv_dom_ll_1_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_srv_dom_ll_f[,1], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel_srv_dom_ll_1_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Male")
lines(data$ages, mle_report$sel_srv_dom_ll_m[,1], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

png(filename = file.path(fig_path, "sel_srv_dom_LL_2_age.png"), width = 7, height = 5, units ="in", res = 250)
par(mfrow = c(1,2))
plot(data$ages, sel_srv_dom_ll_10_f, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Female")
lines(data$ages, mle_report$sel_srv_dom_ll_f[,2], lwd = 3, col = "black", lty = 2)
plot(data$ages, sel_srv_dom_ll_10_m, lwd = 3, col = "red", type = "l", xlab = "Age", ylab = "Selectivity", main = "Male")
lines(data$ages, mle_report$sel_srv_dom_ll_m[,2], lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()


## plot Predicted catch
png(filename = file.path(fig_path, "pred_Catch.png"), width = 8, height = 8, units ="in", res = 250)
par(mfrow = c(2,1))
plot(admb_years, catch_fish1, lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "Catch", main = "Longline")
lines(data$years, mle_report$annual_fixed_catch_pred, lwd = 3, col = "black", lty = 2)
plot(admb_years, catch_fish3, lwd = 3, col = "red", type = "l", xlab = "Year", ylab = "Catch", main = "Trawl")
lines(data$years, mle_report$annual_trwl_catch_pred, lwd = 3, col = "black", lty = 2)
legend('bottomright', legend = c("CA", "SSA"), col= c("red","black"), lty = c(1,2), lwd =2)
dev.off()

################################################################
## Observations
################################################################

###
## Survey DOM LL biomass
###
admb_srv_ll_dom <- as.numeric(unlist(strsplit(tem_rep[(grep("Survey Biomass 1",tem_rep)+2)],split=" ")))[-c(1:4)] #Get numbers at age
obs_srv_ll_dom <- as.numeric(unlist(strsplit(tem_rep[(grep("Survey Biomass 1",tem_rep)+3)],split=" ")))[-c(1:4)] #Get numbers at age

png(filename = file.path(fig_path, "compare_srv1_dom_ll.png"), width = 6, height = 6, units ="in", res = 250)
plot(sab_inputs$srv_dom_ll_bio_yrs, obs_srv_ll_dom, lwd = 3, col = "red", type = "p", xlab = "Year", ylab = "Biomass", ylim = c(0,4000))
lines(admb_years, admb_srv_ll_dom, lwd = 3, lty = 2, col = "blue")
lines(data$years[mle_report$srv_dom_ll_bio_indicator == 1], mle_report$pred_srv_dom_ll_bio[mle_report$srv_dom_ll_bio_indicator == 1], lwd = 3, col = "black", lty = 3)
legend('topleft', legend = c("obs","CA", "SSA"), col= c("red","blue","black"), lty = c(1,2), lwd = 2)
dev.off()


