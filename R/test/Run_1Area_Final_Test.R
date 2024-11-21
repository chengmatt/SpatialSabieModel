# Purpose: Run 1-Area model as a final model for comparison 1960
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 9/4/24

# Set up ------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
library(here)
source(here("R", "Utility_Fxns.R"))
source(here("R", "test", "Francis_Reweight.R"))
scenario = "1960_start" # Model scenario
fig_path = here("Figs", "1-Area-1960")
out_path = here("Output", "Final Models", "1-Area-1960")

# Load in data and parameters
data = readRDS(file = here(out_path, "data.RDS"))
parameters = readRDS(file = here(out_path, "parameters.RDS"))
multiple_shoot = TRUE # whether to do multiple shooting
shoot_iter = 7 # number of iterations to multiple shoot for

# Run Models ---------------------------------------------------------------
bins <- list(
  rep(2:31,2),
  rep(seq(41,99,2),2),
  rep(2:31,2),
  rep(2:31,2),
  rep(seq(41,99,2),2)
)

### Base --------------------------------------------------------------------
# Parameterization
# Start year 1960
# 1 region, 2 fishery fleets, and 3 survey fleets
# Mean recruitment
# Estimate all recruitment deviations except for the last year (does not sum to zero)
# Time block for fixed gear fleet selectivity from 1960 - 2015, 2016 - terminal (logistic)
# Time invariant trawl fishery selectivity (gamma)
# Selectivity for fixed gear deltas and trawl are shared across sexes
# No tag data used
# Movement is fixed at identity

# Set up path for outputting model
model_path = here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final")
dir.create(model_path)

srv_sel_first_param_shared_by_sex = F 
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # no spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "off" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = F # fixed movement
est_prop_male_recruit = "off" # fixed sex-ratio


# Map parameters off
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars
parameters$ln_fixed_sel_pars[] = log(1) # get it stuck out of local minima
parameters$ln_trwl_sel_pars[] = log(5) # get it stuck out of local minima
parameters$trans_srv_q[] = log(7e3)

iter_rw <- 1
wgts_mat <- matrix(data = 0, nrow = iter_rw, ncol = 5)
b0_vec <- vector()

for(iter in 1:iter_rw) {
  
  data = readRDS(file = here(out_path, "data.RDS"))
  
  if(iter > 1) {
    data$obs_fixed_catchatage[] = data$obs_fixed_catchatage[] * wgts[1]
    data$obs_fixed_catchatlgth[] = data$obs_fixed_catchatlgth[] * wgts[2]
    data$obs_srv_catchatage[,,,1] = data$obs_srv_catchatage[,,,1] * wgts[3]
    data$obs_srv_catchatage[,,,2] = data$obs_srv_catchatage[,,,2] * wgts[4]
    data$obs_trwl_catchatlgth[] = data$obs_trwl_catchatlgth[] * wgts[5]
  }
  
  mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, 
                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
  mle_spatial$convergence # check convergence, 0 = converged
  
  # Run more newton steps
  try_improve = tryCatch(expr =
                           for(i in 1:3) {
                             g = as.numeric(mle_obj$gr(mle_spatial$par))
                             h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                             mle_spatial$par = mle_spatial$par - solve(h,g)
                             mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                           }, error = function(e){e})
  
  # check post optimization
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
  wgts <- francis_rwgt(mle_report, "single", bin_list = bins)
  
  wgts_mat[iter, ] <- wgts
  b0_vec[iter] <- mle_report$Bzero
  
  # par(mfrow = c(3,2))
  # plot(b0_vec, type = 'l')
  # plot(wgts_mat[,1], type = 'l')
  # plot(wgts_mat[,2], type = 'l')
  # plot(wgts_mat[,3], type = 'l')
  # plot(wgts_mat[,4], type = 'l')
  # plot(wgts_mat[,5], type = 'l')
  # 
  # par(mfrow = c(1,1))
  # plot(mle_report$SSB_yr)
  
}

lines(mle_report$SSB_yr)

sd_report = sdreport(mle_obj)


### Tagging --------------------------------------------------------------------
# Parameterization
# Start year 1960
# 1 region, 2 fishery fleets, and 3 survey fleets
# Mean recruitment
# Estimate all recruitment deviations except for the last year (does not sum to zero)
# Time block for fixed gear fleet selectivity from 1960 - 2015, 2016 - terminal (logistic)
# Time invariant trawl fishery selectivity (gamma)
# Selectivity for fixed gear deltas and trawl are shared across sexes
# No tag data used
# Movement is fixed at identity

# Set up path for outputting model
model_path = here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final")
dir.create(model_path)

srv_sel_first_param_shared_by_sex = F 
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # no spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = c(1978, 1995, 2017) # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = F # fixed movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1 # negative binomial likelihood
data$evaluate_tag_likelihood = 1 # evaluate likelihood

# Map parameters off
map_fixed_pars = set_up_parameters(data = data, parameters = parameters,
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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars
parameters$ln_fixed_sel_pars[] = log(1) # get it stuck out of local minima
parameters$ln_trwl_sel_pars[] = log(5) # get it stuck out of local minima
parameters$trans_srv_q[] = log(7e3)

iter_rw <- 7
wgts_mat <- matrix(data = 0, nrow = iter_rw, ncol = 5)
b0_vec <- vector()

for(iter in 1:iter_rw) {
  
  data = readRDS(file = here(out_path, "data.RDS"))
  data$tag_likelihood = 1 # negative binomial likelihood
  data$evaluate_tag_likelihood = 1 # evaluate likelihood
  
  if(iter > 1) {
    data$obs_fixed_catchatage[] = data$obs_fixed_catchatage[] * wgts[1]
    data$obs_fixed_catchatlgth[] = data$obs_fixed_catchatlgth[] * wgts[2]
    data$obs_srv_catchatage[,,,1] = data$obs_srv_catchatage[,,,1] * wgts[3]
    data$obs_srv_catchatage[,,,2] = data$obs_srv_catchatage[,,,2] * wgts[4]
    data$obs_trwl_catchatlgth[] = data$obs_trwl_catchatlgth[] * wgts[5]
  }
  
  mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, 
                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
  mle_spatial$convergence # check convergence, 0 = converged
  
  # Run more newton steps
  try_improve = tryCatch(expr =
                           for(i in 1:3) {
                             g = as.numeric(mle_obj$gr(mle_spatial$par))
                             h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                             mle_spatial$par = mle_spatial$par - solve(h,g)
                             mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                           }, error = function(e){e})
  
  # check post optimization
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
  wgts <- francis_rwgt(mle_report, "single", bin_list = bins)
  
  wgts_mat[iter, ] <- wgts
  b0_vec[iter] <- mle_report$Bzero
  
  par(mfrow = c(3,2))
  plot(b0_vec, type = 'l')
  plot(wgts_mat[,1], type = 'l')
  plot(wgts_mat[,2], type = 'l')
  plot(wgts_mat[,3], type = 'l')
  plot(wgts_mat[,4], type = 'l')
  plot(wgts_mat[,5], type = 'l')
  
}

plot(mle_report$SSB_yr)
plot(mle_report$recruitment_yr)

sd_report = sdreport(mle_obj)

