# Purpose: Run 5-Area model as a final model for comparison 1960
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

scenario = "1960_start" # Model scenario
fig_path = here("Figs", "5-Area-1960")
out_path = here("Output", "Final Models", "5-Area-1960")

multiple_shoot = TRUE # whether to do multiple shooting
shoot_iter = 10 # number of iterations to multiple shoot for

# Run Models ---------------------------------------------------------------
# Base Parameterization
# Share survey and fishery selex delta parameters by sex
# Fishery and survey selex is time-invariant
# Don't estimate last 2 years of recruitment
# Tag likelihood is Poisson
# Movement is time-and age-invariant
# Catchability is constant

# The following models build on:
# 1) how movement is parameterized - i.e., age vs. time vs. agextime varying in blocks
# 2) tag likelihoods (poisson vs. negative binomial vs. multinomial)
# 3) Reporting rates (decadal, block via fishery, or spatially reporting)
# and 4) whether spatial q is used for abundane indices

### No Time Blocks -----------------------------------------------------------

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Base")
dir.create(model_path)

# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))

# Set up time-varying movement (data inputs)
data$n_movement_age_blocks = 1 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(0, length(17:31))) 

data$n_movement_time_blocks = 1 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(0, length(1991:2000)),
                                       rep(0, length(2001:2010)), rep(0, length(2011:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src"))     dyn.load(dynlib('TagIntegrated'))          mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### Poisson (Movement) -----------------------------------------------------------------
# Parameterization
# Start year 1960
# 5 regions, 2 fishery fleets, and 2 survey fleets
# Mean recruitment regional
# Estimate all regional recruitment deviations except for the last year (does not sum to zero)
# Time block for fixed gear fleet selectivity from 1960 - 2015, 2016 - terminal (logistic)
# Time invariant trawl fishery selectivity (gamma)
# Selectivity for fixed gear deltas and trawl are shared across sexes
# Tag data used, Poisson
# Movement is estimated (time-varying)
# Spatially-constant q
# Constant reporting rates

#### 4 Time Blocks -----------------------------------------------------------

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Time_4_Move")
dir.create(model_path)

# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))

# Set up time-varying movement (data inputs)
data$n_movement_age_blocks = 1 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(0, length(17:31))) 

data$n_movement_time_blocks = 4 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2000)),
                                       rep(2, length(2001:2010)), rep(3, length(2011:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

#### 3 Time Blocks -----------------------------------------------------------

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Time_3_Move")
dir.create(model_path)

# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))

# Set up time-varying movement (data inputs)
data$n_movement_age_blocks = 1 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(0, length(17:31))) 

data$n_movement_time_blocks = 3 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1995)), rep(1, length(1996:2010)), rep(2, length(2011:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))


#### 2 Time Blocks -----------------------------------------------------------

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Time_2_Move")
dir.create(model_path)
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))

# Set up time-varying movement (data inputs)
data$n_movement_age_blocks = 1 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(0, length(17:31))) 

data$n_movement_time_blocks = 2 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 2 Age Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_2_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 2 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(1, length(17:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))

    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 3 Age Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_3_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 2 Age Blocks x 4 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_2_Time_4_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 2 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(1, length(17:31))) 

data$n_movement_time_blocks = 4 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2000)),
                                       rep(2, length(2001:2010)), rep(3, length(2011:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 3 Age Blocks x 4 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_3_Time_4_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

data$n_movement_time_blocks = 4 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2000)),
                                       rep(2, length(2001:2010)), rep(3, length(2011:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 2 Age Blocks x 3 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_2_Time_3_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 2 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(1, length(17:31))) 

data$n_movement_time_blocks = 3 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1995)), rep(1, length(1996:2010)), rep(2, length(2010:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 3 Age Blocks x 3 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_3_Time_3_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

data$n_movement_time_blocks = 3 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1995)), rep(1, length(1996:2010)), rep(2, length(2010:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 2 Age Blocks x 2 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_2_Time_2_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 2 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:16)), rep(1, length(17:31))) 

data$n_movement_time_blocks = 2 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

### 3 Age Blocks x 2 Time Blocks -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-01-Age_3_Time_2_Move")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

data$n_movement_time_blocks = 2 # time blocks
data$movement_time_block_indicator = c(rep(0, length(1960:1990)), rep(1, length(1991:2021)))

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## Negative Binomial -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-02-NegBin")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "constant" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    compile("TagIntegrated.cpp")
    dyn.load(dynlib('TagIntegrated'))

    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## Multinomial -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-02-Mltnml")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

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
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 2

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Map parameters off
na_map = fix_pars(par_list = parameters, pars_to_exclude = "ln_tag_phi") # fix negative binomial parameter
map_fixed_pars = set_up_parameters(data = data,
                                   parameters = parameters,
                                   na_map = na_map,
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

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## 3 Age Blocks (Reporting - Decadal) -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-03-Decadal")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = c(1978, 1990, 2000, 2010) # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## 3 Age Blocks (Reporting - Fish Block) -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-03-FishBlock")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = c(1978, 1995, 2017) # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## 3 Age Blocks (Reporting - Spatial) -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-03-Space")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = F # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = "space" # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## 3 Age Blocks (Reporting - Fish Block, Spatial Q) -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-04-SptQ")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 3 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(1, length(8:15)), rep(2, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

srv_sel_first_param_shared_by_sex = F
srv_sel_second_param_shared_by_sex = T
fixed_sel_first_param_shared_by_sex  = F
fixed_sel_second_param_shared_by_sex   = T
trwl_sel_first_param_shared_by_sex  = F
trwl_sel_second_param_shared_by_sex  = T
recruit_dev_years_not_to_estimate = 2020:2021 # don't estimate rec devs for last 2 years
srv_q_spatial = T # spatial q
est_init_F = F # estimate initial F
tag_reporting_rate = c(1978, 1995, 2017) # no tag reporting
est_catch_sd = F # fixed catch sd
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# # Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
as.vector(parameters$ln_srv_sel_pars)
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))

## Time and Age-invariant Movement (No Tag Data, No Reporting Spatial Q) -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

# Set up path for outputting model
model_path = here("Output", "Final Models", "5-Area-1960", "5-Area-1960-00-NoTag")
dir.create(model_path)

# Set up age-varying movement (data inputs)
data$n_movement_age_blocks = 1 # two age blocks
data$movement_age_block_indicator = c(rep(0, length(2:7)), rep(0, length(8:15)), rep(0, length(16:31))) 

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))
# set diagnoals equal to 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/") # renormalise

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      parameters$transformed_movement_pars[,i,t,a] = simplex(move_matrix[i,])
    } # end a block
  } # end t block
} # end i region

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
est_movement = T # estimated movement
est_prop_male_recruit = "off" # fixed sex-ratio
data$tag_likelihood = 0
data$tag_likelihood = 1

# Fixed Gear Fishery Selectivity
# data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer") # logistic selectivity for males and females
# data$fixed_sel_by_year_indicator = as.vector(rep(0, length(data$years)), mode = "integer")
# parameters$ln_fixed_sel_pars = array(0, dim = c(length(unique(data$fixed_sel_by_year_indicator)), 2, 2))
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

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

# Share deltas by fishery (estimate sex-specific but share across blocks)
ln_fixed_sel_pars = factor(c(1,2,3,3,4,5,3,3))
# # Share deltas by survey (estimate sex-specific but share across fleets)
ln_srv_sel_pars = factor(c(1,2,3,4,5,2,6,4))
as.vector(parameters$ln_srv_sel_pars)
map_fixed_pars$ln_fixed_sel_pars = ln_fixed_sel_pars
map_fixed_pars$ln_srv_sel_pars = ln_srv_sel_pars
parameters$ln_fixed_sel_pars[] = log(1)
parameters$ln_trwl_sel_pars[] = log(5)

if(multiple_shoot == TRUE) {
  for(i in 1:shoot_iter) {
    if(i > 1) {
      parameters$trans_rec_dev[] = mle_report$recruitment_devs
      parameters$ln_init_rec_dev[] = log(mle_report$init_rec_dev)
      parameters$ln_tag_phi = log(mle_report$tag_phi)
      parameters$ln_trwl_sel_pars[] = log(mle_report$trwl_sel_pars)
      parameters$ln_srv_sel_pars[] = log(mle_report$srv_sel_pars)
      parameters$trans_srv_q[] = log(mle_report$srv_q)
      parameters$logistic_tag_reporting_rate[] = logit(mle_report$tag_reporting_rate)
    }
    
    # Load in tagintegrated model
    setwd(here("src"))
    dyn.load(dynlib('TagIntegrated'))
    
    mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
    pre_optim_sanity_checks(mle_obj)
    
    mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                         gradient  = mle_obj$gr,
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
    post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_obj$env$last.par.best)
    mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
    if(max(abs(mle_obj$gr(mle_obj$env$last.par.best))) < 0.001) break
    
  } # end i
} else {
  
  setwd(here("src")); dyn.load(dynlib('TagIntegrated')); mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T)
  pre_optim_sanity_checks(mle_obj)
  
  mle_spatial = nlminb(start = mle_obj$par, objective = mle_obj$fn,
                       gradient  = mle_obj$gr,
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
  post_optim_sanity_checks(mle_obj = mle_obj, mle_pars = mle_spatial$par)
  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
}

sd_report = sdreport(mle_obj)

# Save stuff
saveRDS(data, file.path(model_path, "data.RDS"))
saveRDS(parameters, file.path(model_path, "parameters.RDS"))
saveRDS(mle_report, file.path(model_path, "mle_report.RDS"))
saveRDS(sd_report, file.path(model_path, "sd_report.RDS"))
saveRDS(mle_spatial, file.path(model_path, "mle_optim.RDS"))
saveRDS(map_fixed_pars, file.path(model_path, "map_fixed_pars.RDS"))
mle_param_list = mle_obj$env$parList(par = mle_obj$env$last.par.best)
saveRDS(mle_param_list, file.path(model_path, "mle_par_list.RDS"))