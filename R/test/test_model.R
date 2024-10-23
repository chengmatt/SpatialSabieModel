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

# Read in data
data = readRDS(here(out_path, "data_5area_agesexmove.RDS"))
parameters = readRDS(here(out_path, "parameters_5area_agesexmove.RDS"))



### Base --------------------------------------------------------------------
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
data$tag_likelihood = 0
data$age_based_movement = 100 # no age based movement
data$sex_based_movement = 1 # sex based movement


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
# mod = readRDS(here(out_path, "5-Area-1960-01-Poisson", "mle_report.RDS"))

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