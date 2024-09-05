# Purpose: Run 5-Area model as a final model for comparison 1977
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
scenario = "1977_start" # Model scenario
fig_path = here("Figs", "5-Area-1977")
out_path = here("Output", "Final Models", "5-Area-1977")

# Load in data and parameters
data = readRDS(file = here(out_path, "data.RDS"))
parameters = readRDS(file = here(out_path, "parameters.RDS"))

# Run Models ---------------------------------------------------------------

### Base --------------------------------------------------------------------
# Parameterization
# Start year 1977
# 5 regions, 2 fishery fleets, and 2 survey fleets
# Mean recruitment regional
# Estimate all regional recruitment deviations except for the last year (does not sum to zero)
# Time block for fixed gear fleet selectivity from 1977 - 2015, 2016 - terminal (logistic)
# Time invariant trawl fishery selectivity (gamma)
# Selectivity for fixed gear deltas and trawl are shared across sexes
# Tag data used, Neg Bin
# Movement is estimated 
# Spatially varying q
# Time block reporting rates

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1977", "5-Area-1977-Base")
dir.create(model_path)

srv_sel_first_param_shared_by_sex = F 
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = T # no spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = c(1978, 1990, 2000, 2010) # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Map parameters off
# na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
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

mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="SpatialSablefishAssessment_TMBExports", hessian = T)
pre_optim_sanity_checks(mle_obj)

mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn, gradient  = mle_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
mle_spatial$convergence # check convergence, 0 = converged

# Run more newton steps
try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(mle_obj$gr(mle_spatial$par))
                           h = optimHess(mle_spatial$par, fn = mle_obj$fn, gr = mle_obj$gr)
                           mle_spatial$par = mle_spatial$par - solve(h,g)
                           mle_spatial$objective = mle_obj$fn(mle_spatial$par)
                         }, error = function(e){e})

# check post optimization
post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
saveRDS(region_key, file.path(model_path, "region_key.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

