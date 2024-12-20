# Purpose: Run sensitivity tests for 5-area spatial models
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 12/17/24

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

multiple_shoot = FALSE # whether to do multiple shooting
shoot_iter = 10 # number of iterations to multiple shoot for

## 3 Age Blocks (Reporting - Fish Block) Sensitivity on ISS -----------------------------------------------------------
# Load in data and parameters
data = readRDS(file = here(out_path, "data_5area_agemove.RDS"))
parameters = readRDS(file = here(out_path, "parameters_5area_agemove.RDS"))

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
parameters$ln_srv_sel_pars[] = log(0.5) # starting this elswhere so it doesn't got off the rails

# Vary ISS
iss <- c(seq(10, 30, 10), seq(50, 200, 10))

for(iter in 1:length(iss)) {
  
  # Set up path for outputting model
  model_path <- here("Output", "Final Models", "5-Area-1960", paste("5-Area-1960-03-FishBlock", iss[iter], sep = "_"))
  dir.create(model_path)
  
  # Vary input ISS here
  for(r in 1:data$n_regions) {
    # Get fixed gear samples
    fxed_gr_ages <- data$obs_fixed_catchatage[,r,which(data$fixed_catchatage_indicator[r,] == 1)] 
    data$obs_fixed_catchatage[,r,which(data$fixed_catchatage_indicator[r,] == 1)]  <- (fxed_gr_ages / colSums(fxed_gr_ages)) * iss[iter] # input varied iss in
    # Get fixed gear samples len
    fxed_gr_lens <- data$obs_fixed_catchatlgth[,r,which(data$fixed_catchatlgth_indicator[r,] == 1)] 
    data$obs_fixed_catchatlgth[,r,which(data$fixed_catchatlgth_indicator[r,] == 1)]  <- (fxed_gr_lens / colSums(fxed_gr_lens)) * iss[iter] # input varied iss in
    # Get trawl gear samples
    trwl_gr_lens <- data$obs_trwl_catchatlgth[,r,which(data$trwl_catchatlgth_indicator[r,] == 1)] 
    data$obs_trwl_catchatlgth[,r,which(data$trwl_catchatlgth_indicator[r,] == 1)]  <- (trwl_gr_lens / colSums(trwl_gr_lens) * iss[iter]) # input varied iss in
    # Get jp survey gear samples
    jp_gr_ages <- data$obs_srv_catchatage[,r,which(data$srv_catchatage_indicator[r,,1] == 1),1]
    data$obs_srv_catchatage[,r,which(data$srv_catchatage_indicator[r,,1] == 1),1] <- (jp_gr_ages / colSums(jp_gr_ages) * iss[iter]) # input varied iss in
    # Get us survey gear samples
    us_gr_ages <- data$obs_srv_catchatage[,r,which(data$srv_catchatage_indicator[r,,2] == 1),2]
    data$obs_srv_catchatage[,r,which(data$srv_catchatage_indicator[r,,2] == 1),2] <- (us_gr_ages / colSums(us_gr_ages) * iss[iter]) # input varied iss in
  } # end r loop
  
  # Map parameters off
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

  # Load in tagintegrated model
  setwd(here("src"))
  dyn.load(dynlib('TagIntegrated'))
  
  mle_obj = MakeADFun(data, parameters, map = map_fixed_pars, DLL="TagIntegrated", hessian = T, silent = T)
  
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

  mle_report = mle_obj$report(mle_obj$env$last.par.best) # get report
  
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
  
  print(iter)
} # end iteration for iss loop
