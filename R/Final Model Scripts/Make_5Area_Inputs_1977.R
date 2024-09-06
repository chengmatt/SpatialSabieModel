# Purpose: To make a universal script that creates data inputs for 5 area models with a model
# start year of 1977 Taken largely from Craig Marsh's work.
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 9/4/24

# Set up ------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
library(here)
source(here("R", "Utility_Fxns.R"))
scenario = "1977_start" # Model scenario

dat_path = here("Data", "5-Area")
fig_path = here("Figs", "5-Area-1977")
out_path = here("Output", "Final Models", "5-Area-1977")
if(!dir.exists(fig_path)) dir.create(fig_path)
if(!dir.exists(out_path)) dir.create(out_path)

### ADMB Data  --------------------------------------------------------------
sab_curr = dget(here("Data", "ADMB", "tem.rdat")) 
sab_ctl = digest_sab_ctl(here("Data", "ADMB", "tem.ctl")) 
sab_inputs = digest_sab_input(here("Data", "ADMB", "tem_2022_na_wh.dat")) 
sab_par = readLines(here("Data", "ADMB", "tem.par")) 
sab_rep = readLines(here("Data", "ADMB", "sable.rep")) 
tem_rep = readLines(here("Data", "ADMB", "tem.rep")) 
admb_out = readLines(here("Data", "ADMB", "SABLE_SARA.dat"))

### Abundance Indices -------------------------------------------------------
design_survey_index = readRDS(file = file.path("Data", "Survey", "regional_abundance_estimates.RDS"))
design_survey_index$area_lab = design_survey_index$NPFMC.Sablefish.Management.Area
design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "Aleutians" ~ "AI",
                                                                 area_lab == "Bering Sea" ~ "BS",
                                                                 area_lab == "Western Gulf of Alaska" ~ "WGOA",
                                                                 area_lab == "Central Gulf of Alaska" ~ "CGOA",
                                                                 area_lab == "Eastern Gulf of Alaska" ~ "EGOA",
                                                                 TRUE ~ area_lab))
design_survey_index = design_survey_index %>% group_by(Country, Year, area_lab) %>% summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)

### Composition Data --------------------------------------------------------
AF_direct_ageing = T # Whether direct ageing is used
max_N_eff = 500/5 # Maximum input sample size
min_N_eff = 50/5  # Minimum input sample size
N_eff_multiplier = 0.4 # Multiplier for input sample size

# read in composition data
fixed_gear_AF_direct = readRDS(file = here(dat_path, "Observer_direct_AF_w_eff.RDS"))
trawl_observer_LF = readRDS(file = here(dat_path, "Observer_trawl_LF_w_eff.RDS"))
fixed_observer_LF = readRDS(file = here(dat_path, "Observer_fixed_LF_w_eff.RDS"))
survery_AF_direct = readRDS(file = here(dat_path, "Survey_direct_AF_w_eff.RDS"))  %>% mutate(country = ifelse(year <= 1993, "Japan", "United States"))

# Make sure input sample sizes are not overly large or small (cap at 500, floor at 80, and multiply by 0.4 if not larger than the max cap)
fixed_gear_AF_direct = fixed_gear_AF_direct %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
survery_AF_direct = survery_AF_direct %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
trawl_observer_LF = trawl_observer_LF %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
fixed_observer_LF = fixed_observer_LF %>% mutate(eff_N = ifelse(eff_N * N_eff_multiplier > max_N_eff, max_N_eff * N_eff_multiplier, eff_N * N_eff_multiplier))
fixed_gear_AF_direct = fixed_gear_AF_direct %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
survery_AF_direct = survery_AF_direct %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
trawl_observer_LF = trawl_observer_LF %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))
fixed_observer_LF = fixed_observer_LF %>% mutate(eff_N = ifelse(eff_N <= min_N_eff, min_N_eff, eff_N))

### Tag Data ----------------------------------------------------------------
tag_recovery_df = readRDS(file = here(dat_path, "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = here(dat_path, "Tag_release_summarised.RDS"))

### Catch Data --------------------------------------------------------------
full_catch_df = readRDS(file = here(dat_path, "Catch_by_year_area_gear.RDS"))
fixed_gear_with_imputation = readRDS(file = here("Data", "5-Area", "fixed_gear_with_imputations_S1.RDS"))
trawl_gear_with_imputation = readRDS(file = here("Data", "5-Area", "trawl_gear_with_imputations_S1.RDS"))

# Model Dimensions --------------------------------------------------------
years = 1977:2021 # Years for spatial model
yrs_to_add = 1977 - 1960 # years for catch imputation
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))
M = 0.104884 # natural mortality rate
reg_lvls = region_key$area[region_key$TMB_ndx + 1] ## need to be in correct order
year_lvls = years ## need to be in correct order

# age bins (males then females)
min_age = 2 # minimum age
max_age = 31 # maximum age
ages = min_age:max_age # age range
sex_age_lvls = paste0(rep(c("M","F"), each = length(ages)), ages)
fixed_gear_AF_direct$sex_age = factor(fixed_gear_AF_direct$sex_age, levels = sex_age_lvls, ordered = T) # create sex length names
survery_AF_direct$sex_age = factor(survery_AF_direct$sex_age, levels = sex_age_lvls, ordered = T) # create sex length names

# length bins (males then females)
length_bins = seq(from = 41, to = 99,by = 2) # number of bins
sex_length_lvls = paste0(rep(c("M","F"), each = length(length_bins)), length_bins)
trawl_observer_LF$sex_length = factor(trawl_observer_LF$sex_length, levels = sex_length_lvls, ordered = T) # create sex length names 
fixed_observer_LF$sex_length = factor(fixed_observer_LF$sex_length, levels = sex_length_lvls, ordered = T) # create sex length names

# TMB Data ----------------------------------------------------------------
data <- list()
data$model = "TagIntegrated" # tag integrated model
data$ages = ages # number of ages
data$years = years # number of years
data$length_bins = length_bins # length bins
data$n_regions = nrow(region_key) # regions
survey_labels = c("Japanese", "Domestic")
data$n_surveys = length(survey_labels) # number of surveys
data$n_projections_years = 0 # projeciton years
data$do_projection = 0 # don't do projections

# Dimensions 
n_regions = data$n_regions
n_ages = length(data$ages) 
n_length_bins = length(data$length_bins) # the last length bin value is the minimum for a length plus group
n_projyears = length(data$years) +  data$n_projections_years
n_years = length(data$years)
projyears = min(data$years):(max(data$years) + data$n_projections_years)


# Projections
data$future_recruitment_type = 0 # simulate from lognormal
data$year_ndx_for_empirical_resampling = c(0,n_years - 1) # where to resmaple from (years)
data$future_fishing_type = 0 # user supplifed F
data$future_fishing_inputs_trwl = array(0.1, dim = c(data$n_regions, data$n_projections_years)) # specified F
data$future_fishing_inputs_fixed = array(0.1, dim = c(data$n_regions, data$n_projections_years)) # specific F

### Recruitment Specifications ----------------------------------------------
data$global_rec_devs = 0 # regional rec devs
data$rec_devs_sum_to_zero = 0 # rec devs do not sum to zero
data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(n_years) # vector for summing to zero (not used)
data$n_init_rec_devs = sab_inputs$yrs[1] - sab_ctl$rec_styr - 2 # number of non-equilibrium age-structure
data$spawning_time_proportion = rep(0, n_projyears) # when is spawning caclulated (start of year)
data$sigma_R = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("sigr:",sab_par)+1],split=" ")))) # sigma R used
data$SrType = 3 # mean recruitment

### Biology Specifications --------------------------------------------------
data$ageing_error_matrix = sab_inputs$ageing_error # ageging error matrix
data$M = matrix(0.104884, nrow = n_ages, ncol = n_projyears) # natural mortality
data$maturity = matrix(sab_curr$growthmat[,5], nrow = n_ages, ncol = n_projyears, byrow = F) # maturity

# Weight at age 
weight_at_age_f = as.numeric(unlist(strsplit(sab_rep[(grep("Weight Females",sab_rep)+1):(grep("Weight Females",sab_rep)+n_years)],split=" "))) # female weight at age
weight_at_age_m = as.numeric(unlist(strsplit(sab_rep[(grep("Weight Males",sab_rep)+1):(grep("Weight Males",sab_rep)+n_years)],split=" "))) # male weight at age
weight_at_age_f = subset(weight_at_age_f, subset = !is.na(weight_at_age_f))
weight_at_age_m = subset(weight_at_age_m, subset = !is.na(weight_at_age_m))
weight_at_age_f_mat = matrix(weight_at_age_f, byrow = T, nrow = n_years, ncol = n_ages) # turn to matrix
weight_at_age_m_mat = matrix(weight_at_age_m, byrow = T, nrow = n_years, ncol = n_ages) # turn to matrix
data$female_mean_weight_by_age = t(weight_at_age_f_mat)
data$male_mean_weight_by_age = t(weight_at_age_m_mat)

# Age length transitions
data$female_age_length_transition = data$male_age_length_transition = array(0, dim = c(n_ages, n_length_bins, n_years), dimnames = list(data$ages, data$length_bins, data$years))
data$male_age_length_transition[,,as.character(min(years):1994)] = sab_inputs$male_length_at_age_60_96
data$female_age_length_transition[,,as.character(min(years):1994)] = sab_inputs$female_length_at_age_60_96
data$male_age_length_transition[,,as.character(1995:max(years))] = sab_inputs$male_length_at_age_97_22
data$female_age_length_transition[,,as.character(1995:max(years))] = sab_inputs$female_length_at_age_97_22

### Movement Specifications -------------------------------------------------
data$n_movement_time_blocks = 1 # 1 movement block (fixed)
data$apply_fixed_movement = 0 # estimated movement matrix
data$fixed_movement_matrix = array(0, dim = c(n_regions, n_regions, data$n_movement_time_blocks)) # fixed movement
data$movement_matrix = array(0, dim = c(n_regions, n_regions, data$n_movement_time_blocks)) # estiamted movement matrix
data$movement_time_block_indicator = rep(0, n_years) # time block index for movement
movement_matrix = fixed_movement_matrix = matrix(0, nrow = n_regions, ncol = n_regions)
diag(fixed_movement_matrix) = 1 # set 1 for fixed movement matrix
diag(movement_matrix) = 0.9 # set 0.9 for estimated movement matrix (cant be 1)
data$fixed_movement_matrix[,,1] = fixed_movement_matrix # define movement fixed 1 region
movement_matrix = movement_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1) # random deviate here for 0s in movement matrix
# Normalize so rows sum to 1 in estimated movement matrix
movement_matrix = sweep(movement_matrix, 1, STATS = rowSums(movement_matrix), "/")
data$movement_matrix[,,1] = movement_matrix
data$do_recruits_move = 0 # recruitment is applied after movement occurs
data$tag_likelihood = 0 # Poisson likelihood
data$evaluate_tag_likelihood = 1 # evaluate tag data

### Fishing Specifications --------------------------------------------------
data$F_method = 1 # Use hybrid method (i.e., Fs not estimated)
data$F_max = 3 # max F for hybrid method
data$F_iterations = 4 # Iterations for hybrid method
data$prop_F_hist = 0 # historical F applied to initial age structure

# Munging for catch data
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

### Fishery and Survey Selectivity -----------------------------------------------------
# Fixed Gear Fishery Selectivity
data$fixed_sel_type = as.vector(rep(0, 2), mode = "integer") # logistic selectivity for males and females
# selectivity time blocks (n = 2)
data$fixed_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")
data$fixed_sel_by_year_indicator[projyears > 2015] = 1 # time block 2

# Trawl Gear Fishery Selectivity
data$trwl_sel_type = as.vector(rep(1, 1), mode = "integer") # gamma selectuivity for males and females
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer") # Single time block

# Survey Selectivity
data$srv_sel_type = matrix(0, nrow = 1, ncol = data$n_surveys) # logistic 
data$srv_sel_by_year_indicator = matrix(0, nrow = n_projyears, ncol = data$n_surveys) # blocks for survey

### Tagging -----------------------------------------------------------------
# Tag releases from 1979 - 2016
tag_release_years = c(1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016)
data$tag_release_event_this_year = rep(0, n_years) # initialize vector
data$tag_release_event_this_year[data$years %in% tag_release_years] = 1 # do tag releases in the years specified above
data$male_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))
data$female_tagged_cohorts_by_age = array(0, dim = c(n_ages, n_regions, length(tag_release_years)))
data$n_years_to_retain_tagged_cohorts_for = 15 # retain tagged cohorts for 15 years (dump after) 
data$initial_tag_induced_mortality = rep(0.1, length(tag_release_years)) # tag mortality
data$annual_tag_shedding_rate = 0.02 # tag shedding

# Loop through to fill in stuff
for(y_ndx in 1:length(tag_release_years)) {
  this_release_year = tag_release_df %>% filter(release_year == tag_release_years[y_ndx])
  if(nrow(this_release_year) > 0) {
    regions_to_release = unique(this_release_year$region_release) 
    for(r_ndx in 1:length(regions_to_release)) {
      region_ndx = region_key$TMB_ndx[which(region_key$area %in% regions_to_release[r_ndx])] + 1
      data$male_tagged_cohorts_by_age[, region_ndx, y_ndx] = (this_release_year %>% filter(sex == "M", region_release == regions_to_release[r_ndx]))$Nage_at_release
      data$female_tagged_cohorts_by_age[, region_ndx, y_ndx] = (this_release_year %>% filter(sex == "F", region_release == regions_to_release[r_ndx]))$Nage_at_release
    } # end r ndx
  } else {
    ## no tag release event
    data$tag_release_event_this_year[data$years %in% tag_release_years[y_ndx]] = 0
  }
} # end y ndx

# Munging tag recovery stuff
include_tag_recoveries = T # tags are turned off using evaluate-tag-likelihood switch
include_zero_tag_recovery_events = T
tag_recovery_years = 1978:2016 # tag recovery years
data$tag_recovery_indicator_by_year = rep(0, n_years) ## no tag releases
data$obs_tag_recovery = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))
data$tag_recovery_indicator = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))

if(include_tag_recoveries) {
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

## are there observations for the plus group
if(sum(data$tag_recovery_indicator[1:5,,]) != 0)
  cat("found recoveries in the first year of release\n")
if(sum(data$tag_recovery_indicator[(dim(data$tag_recovery_indicator)[1] - data$n_regions + 1):dim(data$tag_recovery_indicator)[1],,]) != 0)
  cat("found recoveries in the plus group\n")
cat("The number of recovery events ", sum(data$tag_recovery_indicator), "\n")


### Compositions ------------------------------------------------------------
# fixed gear age frequencies
fixed_gear_AF = fixed_gear_AF_direct %>% filter(year %in% years) # filter to the correct years
fixed_gear_AF$region_f = factor(fixed_gear_AF$region, levels = reg_lvls, ordered = T) # relevel regions
fixed_gear_AF$year_f = factor(fixed_gear_AF$year, levels = year_lvls, ordered = T) # relevel years

# trawl gear length frequencies
trawl_observer_LF = trawl_observer_LF %>% filter(year %in% years) # filter to correct years
trawl_observer_LF$region_f = factor(trawl_observer_LF$region, levels = reg_lvls, ordered = T) # relevel regions
trawl_observer_LF$year_f = factor(trawl_observer_LF$year, levels = year_lvls, ordered = T) # relevel years

# fixed gear length frequencies
fixed_observer_LF = fixed_observer_LF %>% filter(year %in% years) # filter to the correct years
fixed_observer_LF$region_f = factor(fixed_observer_LF$region, levels = reg_lvls, ordered = T) # relevel regions
fixed_observer_LF$year_f = factor(fixed_observer_LF$year, levels = year_lvls, ordered = T) # relelve years

# survey age frequnecies
survery_AF = survery_AF_direct %>% filter(year %in% years) # filter to ccorrect years
survery_AF$region_f = factor(survery_AF$region, levels = reg_lvls, ordered = T) # relelve regions
survery_AF$year_f = factor(survery_AF$year, levels = year_lvls, ordered = T) # relevel years

# create indicator variables
fixed_AF_indicator = fixed_gear_AF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
trawl_LF_indicator = trawl_observer_LF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)
fixed_LF_indicator = fixed_observer_LF %>% filter(!year %in% unique(fixed_gear_AF$year)) %>% ungroup() %>% dplyr::select(region_f, year_f, P) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(P), values_fn = indicator_fun)
survey_AF_indicator = survery_AF %>% ungroup() %>% dplyr::select(region_f, year_f, year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(year), values_fn = indicator_fun)

# fixed gear fishery age frequencies
data$fixed_catchatage_covar_structure = 0 # iid
data$fixed_catchatage_comp_likelihood = 0 # ADMB multinomial 
data$fixed_catchatage_indicator = as.matrix(fixed_AF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years) # indicator variable for fixed gear
data$obs_fixed_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years)) # fill in observations for fixed gear catch
obs_years = unique(fixed_gear_AF$year) # years
obs_reg = unique(fixed_gear_AF$region) # regions
for(y_ndx in 1:length(obs_years)) {
  this_year_ndx = which(years %in% obs_years[y_ndx])
  for(r_ndx in 1:length(obs_reg)) {
    this_region_ndx = region_key$TMB_ndx[which(region_key$area %in% obs_reg[r_ndx])] + 1
    this_df = (fixed_gear_AF %>% filter(region == obs_reg[r_ndx], year == obs_years[y_ndx]))
    if(nrow(this_df) > 0)
      data$obs_fixed_catchatage[,this_region_ndx, this_year_ndx] = this_df$P * this_df$eff_N # propagate numbers by age and sex
  }
}

# trawl gear fishery length frequencies
data$trwl_catchatlgth_covar_structure = 0 # iid
data$trwl_catchatlgth_comp_likelihood = 0 # ADMB multinomial
data$trwl_catchatlgth_indicator = as.matrix(trawl_LF_indicator %>% dplyr::select(-region_f), nrow = n_regions, ncol = n_years) # indicator
data$obs_trwl_catchatlgth = array(0, dim = c(n_length_bins * 2, n_regions, n_years)) # 
obs_years = sort(unique(trawl_observer_LF$year)) # years
obs_reg = unique(trawl_observer_LF$region) # regions
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

# Fixed gear fishery length frequnecies
data$fixed_catchatlgth_covar_structure = 0 # iid
data$fixed_catchatlgth_comp_likelihood = 0 # ADMB multinomial
data$fixed_catchatlgth_indicator = array(0, dim = c(n_regions, n_years))
data$obs_fixed_catchatlgth = array(0, dim = c(n_length_bins * 2, n_regions, n_years))
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

# Survey LL proportions at age
data$srv_catchatage_covar_structure = 0 # iid
data$srv_catchatage_comp_likelihood = rep(0, data$n_surveys) # ADMB multinomial
data$srv_catchatage_indicator = array(0, dim = c(n_regions, n_years,data$n_surveys))
US_survey_years = unique((survery_AF %>% filter(country == "United States"))$year)
jap_survey_years = unique((survery_AF %>% filter(country == "Japan"))$year)
data$obs_srv_catchatage = array(0, dim = c(n_ages * 2, n_regions, n_years, data$n_surveys))
obs_years = unique(survery_AF$year)
obs_reg = unique(survery_AF$region)
countries = unique(survery_AF$country)
# populate container
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

### Abundance Indices -------------------------------------------------------
survey_index = design_survey_index %>% filter(Year %in% years)
survey_index$region_f = factor(survey_index$area_lab, levels = reg_lvls, ordered = T)
survey_index$year_f = factor(survey_index$Year, levels = year_lvls, ordered = T)
survey_index_indicator = survey_index %>% ungroup() %>% dplyr::select(region_f, year_f, Year) %>% complete(year_f, region_f)  %>% pivot_wider(names_from = year_f, values_from = c(Year), values_fn = indicator_fun)

# LL survey abundance indices
data$obs_srv_bio = data$obs_srv_se = array(0, dim = c(n_regions, n_years, data$n_surveys))
data$srv_bio_indicator = array(0, dim = c(n_regions, n_years,data$n_surveys))
countries = unique(design_survey_index$Country)
obs_years = unique(design_survey_index$Year)
obs_reg = unique(design_survey_index$area_lab)
for(c_ndx in 1:length(countries)) {
  for(y_ndx in 1:length(obs_years)) {
    mod_yr_ndx = which(data$years == obs_years[y_ndx])
    for(r_ndx in 1:length(region_key$area)) {
      this_df = design_survey_index %>% filter(Country == countries[c_ndx], Year == obs_years[y_ndx], area_lab == region_key$area[r_ndx])
      if(nrow(this_df) > 0) {
        data$srv_bio_indicator[r_ndx, mod_yr_ndx, c_ndx] = 1
        data$obs_srv_bio[r_ndx, mod_yr_ndx, c_ndx] = this_df$sum_estimates## model numbers are in 1000's is in kilo tonnes
        data$obs_srv_se[r_ndx, mod_yr_ndx, c_ndx] = (this_df$se)  # scale by 1000 as well
      }
    }
  }
}

# likelihoods for abundance indices
data$srv_bio_likelihood = rep(1, data$n_surveys) # lognormal call (some differences in just how se are parameterized, it should be correct like this)
data$srv_obs_is_abundance = rep(1, data$n_surveys) # calculate index by abundance
data$srv_q_by_year_indicator = matrix(0, nrow = n_years, ncol = data$n_surveys)
data$srv_q_transformation = rep(0, data$n_surveys) # log transformation
data$q_is_nuisance = rep(0, data$n_surveys) # estiamted as a free parameter

# TMB Parameters ----------------------------------------------------------
parameters <- list()
parameters$ln_mean_rec = rnorm(data$n_regions, suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_mean_rec:",sab_par)+1],split=" ")))), 0.3)
parameters$ln_sigma_R = log(data$sigma_R) # sigma R
parameters$ln_sigma_init_devs = log(0.4) # initial sigma R
parameters$trans_SR_pars = rep(qlogis(0.8), 1) # SR parameters - not used
parameters$logistic_prop_recruit_male = rep(logit(0.5), length(data$years)) # sex ratio - fixed at 1:1
data$map_simplex_ycs_estimated = rep(0:(length(data$years) - 1))
data$standardise_ycs = 0 # freely estimate

if(data$rec_devs_sum_to_zero == 1) { # if sum to zero
  if(data$global_rec_devs == 1) { # if global rec devs
    parameters$trans_rec_dev = matrix(0, nrow = 1, ncol = n_years - 1)
  } else {
    parameters$trans_rec_dev = matrix(0, nrow = data$n_regions, ncol =  n_years - 1)
  }
} else {
  if(data$global_rec_devs == 1) { # if global rec devs
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = 1, ncol = n_years)
  } else {
    parameters$trans_rec_dev = matrix(0.5 * data$sigma_R^2, nrow = data$n_regions, ncol =  n_years)
  }
}

# Initialize inital age structures
if(data$n_init_rec_devs == 0) {
  parameters$ln_init_rec_dev = 0
} else {
  parameters$ln_init_rec_dev = rep(0, data$n_init_rec_devs)
}

# Selectivity parameters
parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_trwl_sel_pars = array(0, dim = c(length(unique(data$trwl_sel_by_year_indicator )), 2, 2))
parameters$ln_srv_sel_pars = array(0, dim = c(max(apply(data$srv_sel_by_year_indicator, 2,FUN = function(x){length(unique(x))})), 2, 2, data$n_surveys))

# Set up fishery selectivity parameters
parameters$ln_fixed_sel_pars[1,1,1] = parameters$ln_fixed_sel_pars[2,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_m:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,2,1] = parameters$ln_fixed_sel_pars[2,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,1,2] = parameters$ln_fixed_sel_pars[2,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish1_f:",sab_par)+1],split=" "))))
parameters$ln_fixed_sel_pars[1,2,2] = parameters$ln_fixed_sel_pars[2,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish1_f:",sab_par)+1],split=" "))))

# Set up trawl selectivity parameters
parameters$ln_trwl_sel_pars[1,1,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_m:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,1] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,1,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_fish3_f:",sab_par)+1],split=" "))))
parameters$ln_trwl_sel_pars[1,2,2] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_fish3_f:",sab_par)+1],split=" "))))

# Set up survey selectivity parameters
parameters$ln_srv_sel_pars[1,1,1,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_m:",sab_par)+1],split=" "))))
parameters$ln_srv_sel_pars[1,2,1,1:data$n_surveys] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_m:",sab_par)+1],split=" "))))
parameters$ln_srv_sel_pars[1,1,2,1:data$n_surveys] =  suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_a50_srv1_f:",sab_par)+1],split=" "))))
parameters$ln_srv_sel_pars[1,2,2,1:data$n_surveys] = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_delta_srv1_f:",sab_par)+1],split=" "))))

# Set up movement parameters
parameters$transformed_movement_pars = array(NA, dim = c(n_regions - 1, n_regions, data$n_movement_time_blocks))
for(i in 1:n_regions) parameters$transformed_movement_pars[,i,1] = simplex(movement_matrix[i,])
# Set up tagging parameters
parameters$logistic_tag_reporting_rate = matrix(logit(0.50), nrow = data$n_regions, ncol = max(sum(data$tag_recovery_indicator_by_year),1))
parameters$ln_tag_phi = log(0.5) # Negative binomial phi parameter

# Set up q parameters
parameters$trans_srv_q = array(log(7e3), dim =c(n_regions, max(apply(data$srv_q_by_year_indicator, 2,FUN = function(x){length(unique(x))})) , data$n_surveys))

# Fishing mortality and catch parameters
parameters$ln_fixed_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish1:",sab_par)+1],split=" ")))) # fixed gear average
parameters$ln_fixed_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish1:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5)))) # mortality devs
parameters$ln_trwl_F_avg = suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_avg_F_fish3:",sab_par)+1],split=" ")))) # trawl gear average
parameters$ln_trwl_F_devs = log(1.4*exp(aperm(jitter(array(suppressWarnings(as.numeric(unlist(strsplit(sab_par[grep("log_F_devs_fish3:",sab_par)+1],split=" "))))[-1], dim = c(n_projyears, n_regions)), amount = 0.5)))) # mortality devs
parameters$ln_init_F_avg = parameters$ln_fixed_F_avg # initial fixed gear deviations (set equal)
parameters$ln_catch_sd = log(0.02) # deviation of catch

# Data weighitng parameters (Dirichlet-Mutlinomial)
parameters$trans_trwl_catchatlgth_error = log(1)
parameters$trans_fixed_catchatlgth_error = log(1)
parameters$trans_fixed_catchatage_error = log(1)
parameters$trans_srv_catchatage_error = rep(log(1), data$n_surveys)

validate_input_data_and_parameters(data, parameters)

# Save parameters
saveRDS(data, file.path(out_path, "data.RDS"))
saveRDS(parameters, file.path(out_path, "parameters.RDS"))

# Data Input Figures ------------------------------------------------------
# Input observations
ggsave(here(fig_path, "InputObs.png"),
       plot_input_observations(data, region_key = region_key, survey_labels = survey_labels) )

# Time block specifications
ggsave(here(fig_path, "TimeBlock.png"),
       plot_input_timeblocks(data, survey_labels = survey_labels) )

# Catch
ggsave(here(fig_path, "Catch.png"),
       plot_input_catches(data, region_key = region_key) )

# Abundance Indices
ggsave(here(fig_path, "AbdIdx.png"),
       ggplot(design_survey_index, aes(Year, sum_estimates, col = Country, linetype = Country, fill = Country, alpha = 0.3)) +
         geom_ribbon(aes(ymin = LCI, ymax = UCI)) +
         geom_line(linewidth = 1) +
         labs(x = "Year", y = "Relative index (RPN)", linetype = "", col = "", fill = "") +
         guides(alpha = "none") +
         facet_wrap(~area_lab, scales = "free_y") +
         theme_bw())

