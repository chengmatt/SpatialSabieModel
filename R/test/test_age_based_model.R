library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
library(here)
source(here("R", "Utility_Fxns.R"))

scenario = "1960_start" # Model scenario
fig_path = here("Figs", "5-Area-1960")
out_path = here("Output", "Final Models", "5-Area-1960")

# Run Models ---------------------------------------------------------------

### Base --------------------------------------------------------------------
# Read in tag data
tag_recovery_df = readRDS(file = here("Data", "5-Area", "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = here("Data", "5-Area", "Tag_release_summarised.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-age")
dir.create(model_path)

# Load in data and parameters
data = readRDS(file = here(out_path, "data.RDS"))
parameters = readRDS(file = here(out_path, "parameters.RDS"))
multiple_shoot = T # whether to do multiple shooting
shoot_iter = 10 # number of iterations to multiple shoot for
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))

# Parameterization
srv_sel_first_param_shared_by_sex = F 
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "off" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement (time-invariant)
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Other specifications
data$age_based_movement = 1 # age based movement

# Set up tag release years
tag_release_years = c(1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985,
                      1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
                      1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
                      2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
                      2017: 2021)

include_tag_recoveries = T
include_zero_tag_recovery_events = T
tag_recovery_years = 1978:2020

data$tag_recovery_indicator_by_year = rep(0, length(data$years)) ## no tag releases
data$obs_tag_recovery = array(0, dim = c(length(data$ages), data$n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), data$n_regions, length(tag_recovery_years)))
data$tag_recovery_indicator = array(0, dim = c(data$n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), data$n_regions, length(tag_recovery_years)))

if(include_tag_recoveries) {
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
          model_tag_release_ndx = get_tag_release_ndx(release_ndx, release_event_ndx, data$n_regions)
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


# Set up movement parameterization ----------------------------------------
# Time variation
data$n_movement_time_blocks = 1
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(0, length(1990:2021)))

# Age variation
data$n_movement_age_blocks = 3 # three age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
# diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2]) = 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
# renormalise
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/")

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

# Fixed Gear Fishery Selectivity
data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails
parameters$ln_tag_phi = log(2)
validate_input_data_and_parameters(data, parameters)

setwd(here("src"))
compile("TagIntegrated_v1.cpp")
dyn.load(dynlib('TagIntegrated_v1'))
dyn.unload(dynlib('TagIntegrated_v1'))
dyn.load(dynlib('TagIntegrated_v1'))

# Run Model ---------------------------------------------------------------
data$tag_likelihood = 1
# na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_tag_phi")) # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data, 
                                   parameters = parameters,
                                   na_map = NULL,
                                   srv_sel_first_param_shared_by_sex = srv_sel_first_param_shared_by_sex,
                                   srv_sel_second_param_shared_by_sex = srv_sel_second_param_shared_by_sex,
                                   fixed_sel_first_shared_by_sex  = fixed_sel_first_param_shared_by_sex,
                                   fixed_sel_second_shared_by_sex   = fixed_sel_second_param_shared_by_sex,
                                   trwl_sel_first_shared_by_sex  = trwl_sel_first_param_shared_by_sex,
                                   trwl_sel_second_shared_by_sex  = trwl_sel_second_param_shared_by_sex,
                                   recruit_dev_years_not_to_estimate = recruit_dev_years_not_to_estimate,
                                   srv_q_spatial = srv_q_spatial,
                                   tag_reporting_rate = tag_reporting_rate,
                                   est_init_F = est_init_F,
                                   est_catch_sd = est_catch_sd,
                                   est_movement = est_movement,
                                   est_prop_male_recruit = est_prop_male_recruit)

# make ad object
mle_obj = MakeADFun(data, parameters, 
                    map = map_fixed_pars, 
                    DLL="TagIntegrated_v1", hessian = T)

pre_optim_sanity_checks(mle_obj)

mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                     gradient  = mle_obj$gr,
                     control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report

# Run more newton steps
try_improve = tryCatch(expr =
                         for(i in 1:3) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }, error = function(e){e})

# check post optimization
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
sd_report = sdreport(mle_obj)

# Convert TMB code to R (testing) ------------------------------------------
mod = readRDS(here(out_path, "5-Area-1960-01-Poisson", "mle_report.RDS"))

# make movement matrix age specific
nageblock = 3
ages = 1:30
nregion = 5
annual_tag_shedding_rate = 0.02
n_years_to_retain_tagged_cohorts_for = data$n_years_to_retain_tagged_cohorts_for

# make movement matrix
age_block = c(rep(1, 10), rep(2, 10), rep(3, 10))
move_mat_test = mod$fixed_movement_matrix[,,] # from, to, time block, age block
move_mat_test = replicate(nageblock, move_mat_test, simplify = "array")

init_age_f_mat = matrix(0, nrow = length(ages), ncol = nregion)
tagged_natagef = mod$tagged_natage_f
natagef = mod$natage_f
mean_recruit = mod$mean_rec

# need to set up movement matrix
movement_mat = array(1, c(5,5,2,nageblock))
cache_log_k_value = rep(0, 5)
for(t in 1:2) {
  cache_log_k_value = rep(0, 5)
  for(k in 1:nregion) cache_log_k_value[k] = nregion - k
  for(region_ndx in 1:nregion) {
    stick_length = 1.0
    for(ab in 1:3) {
    for(k in 1:nregion) {
        movement_mat[region_ndx,k,t,ab] = stick_length * invlogit(transformed_movement_pars[k, region_ndx, t, ab] - cache_log_k_value[k])
        stick_length = stick_length - movement_mat[region_ndx,k,t,ab]
      } # end k 
      # plus group
      movement_matrix[region_ndx, n_regions - 1, t, ab] = stick_length
    } # end age block
  } # end region
  } # end t

# Initial age structure section
for(age_ndx in 1:length(ages)) {
  init_age_f_mat[age_ndx,] = mod$init_natage_f[age_ndx,] %*% move_mat_test[,,1,1]
}

identical((mod$init_natage_f %*% move_mat_test[,,1,1]), init_age_f_mat)

# Repeat again, for some reason to redo cycle
# Loop through regions and ages
for (age_ndx in 1:length(ages)) { 
  for (to in 1:nregion) {
    # Use sum to replicate matrix multiplication for this age index
    init_age_f_mat[age_ndx, to] <- sum(mod$init_natage_f[age_ndx, ] * move_mat_test[, to, 1, age_block[age_ndx]])
  }
}

reshape2::melt(mle_report$movement_matrix) %>% 
  rename(From = Var1, To = Var2) %>% 
  mutate(
    From = case_when(
      From == 1 ~ "BS",
      From == 2 ~ "AI",
      From == 3 ~ "WGOA",
      From == 4 ~ "CGOA",
      From == 5 ~ "EGOA"
    ), From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), To = factor(To, levels = rev(c("BS", "AI", "WGOA", "CGOA", "EGOA")))) %>% 
  ggplot(aes(x = From, y = To, fill = value, label = round(value, 2))) +
  geom_tile(alpha = 0.55) +
  geom_text() +
  scale_fill_viridis_c() +
  theme_test() +
  facet_grid(Var3~Var4)

plot((mod$init_natage_f %*% move_mat_test[,,1,1])[,3])
lines(init_age_f_mat[,3])

# movement for tagged partition
for(tag_ndx in 1:n_years_to_retain_tagged_cohorts_for) {
  for(release_region_ndx in 1:nregion) {
      tag_release_event_ndx =  tag_ndx * nregion + release_region_ndx # get release event number
      for(age_ndx in 1:length(ages)) {
      tagged_natagef[,,tag_release_event_ndx] = tagged_natagef[age_ndx,,tag_release_event_ndx] %*% move_mat_test[,,1,age_block[age_ndx]]
      tagged_natagef[,,tag_release_event_ndx] = tagged_natagef[,,tag_release_event_ndx] * exp(-annual_tag_shedding_rate)
    } # end age index
  } # end release index
} # end tag index

# Now do movement for the population itself
for(age_ndx in 1:length(ages)) {
  natagef[age_ndx,,2] = natagef[age_ndx,,2] %*% move_mat_test[,,1,age_block[age_ndx]]
}

# if recruits move = 0 (basically just resetting the movement matrix)
for(region_ndx in 1:nregion) {
  natagef[1, region_ndx, 2] = mean_recruit[region_ndx] * devs * prop_recruit
}

# for projection
for(age_ndx in 1:length(ages)) {
  natagef[age_ndx,,proj] = natagef[age_ndx,,proj] %*% move_mat_test[,,1,age_block[age_ndx]]
}

# need to adjust likelihoods to accomodate ages