# Purpose: To conduct francis reweighting 

francis_rwgt <- function(mle_report, area, bin_list) {
  
  if(area == "single") {
    # Get observed quantities
    obs_list <- list(
      obs_fixed_caa = mle_report$obs_fixed_catchatage[,1,mle_report$fixed_catchatage_indicator == 1],
      obs_fixed_cal = mle_report$obs_fixed_catchatlgth[,1,mle_report$fixed_catchatlgth_indicator == 1],
      obs_srvjp_caa = mle_report$obs_srv_catchatage[,1,mle_report$srv_catchatage_indicator[,,1] == 1,1],
      obs_srvdom_caa = mle_report$obs_srv_catchatage[,1,mle_report$srv_catchatage_indicator[,,2] == 1,2],
      obs_srvtrwl_cal = mle_report$obs_trwl_catchatlgth[,1,mle_report$trwl_catchatlgth_indicator == 1]
    )
    
    # Get predicted quantities
    pred_list <- list(
      pred_fixed_caa = mle_report$pred_fixed_catchatage[,1,mle_report$fixed_catchatage_indicator == 1],
      pred_fixed_cal = mle_report$pred_fixed_catchatlgth[,1,mle_report$fixed_catchatlgth_indicator == 1],
      pred_srvjp_caa = mle_report$pred_srv_catchatage[,1,mle_report$srv_catchatage_indicator[,,1] == 1,1],
      pred_srvdom_caa = mle_report$pred_srv_catchatage[,1,mle_report$srv_catchatage_indicator[,,2] == 1,2],
      pred_srvtrwl_cal = mle_report$pred_trwl_catchatlgth[,1,mle_report$trwl_catchatlgth_indicator == 1]
    )
    
    wts <- vector()
    
    for(i in 1:length(obs_list)) {
      
      # Get input sample sizes
      tmp_iss_obs <- apply(obs_list[[i]], 2, sum)
      tmp_obs <- t(apply(X = obs_list[[i]], MARGIN = 2, FUN = function(x) x / sum(x))) # normalize tmp obs 
      tmp_exp <- t(pred_list[[i]]) # get temporrary pred variable
      
      # Set up reweighting vectors
      exp_bar <- vector() # mean expected
      obs_bar <- vector() # mean observed
      v_y <- vector() # variance
      w_denom <- vector() # weight factor in denominator
      
      for(y in 1:dim(tmp_obs)[1]) {
        exp_bar[y] <- sum(bin_list[[i]] * tmp_exp[y,]) # get mean pred comps
        obs_bar[y] <- sum(bin_list[[i]] * tmp_obs[y,]) # get mean obs comps
        v_y[y] <- sum(bin_list[[i]]^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
        w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/200) # get weights (fixing at constant 200)
      } # end y loop
      
      # Get weights 
      wts[i] <- 1 / var(w_denom)
      
    } # end i
  } # if single area
  
  if(area == "spatial") {
    # Get observed quantities
    obs_list <- list(
      obs_fixed_caa = mle_report$obs_fixed_catchatage,
      obs_fixed_cal = mle_report$obs_fixed_catchatlgth,
      obs_srvjp_caa = mle_report$obs_srv_catchatage,
      obs_srvdom_caa = mle_report$obs_srv_catchatage,
      obs_srvtrwl_cal = mle_report$obs_trwl_catchatlgth
    )

    # Get predicted quantities
    pred_list <- list(
      pred_fixed_caa = mle_report$pred_fixed_catchatage,
      pred_fixed_cal = mle_report$pred_fixed_catchatlgth,
      pred_srvjp_caa = mle_report$pred_srv_catchatage,
      pred_srvdom_caa = mle_report$pred_srv_catchatage,
      pred_srvtrwl_cal = mle_report$pred_trwl_catchatlgth
    )
    
    # Get indicators
    idx_list <- list(
      idx_fixed_caa = mle_report$fixed_catchatage_indicator,
      idx_fixed_cal = mle_report$fixed_catchatlgth_indicator,
      idx_srvjp_caa = mle_report$srv_catchatage_indicator[,,1],
      idx_srvdom_caa = mle_report$srv_catchatage_indicator[,,2],
      idx_srvtrwl_cal = mle_report$trwl_catchatlgth_indicator
      )
    
    # set up storage
    wts <- matrix(0, nrow = length(obs_list), ncol = dim(idx_list[[1]])[1])
    
    for(i in 1:length(obs_list)) {
      for(r in 1:dim(idx_list[[i]])[1]) {
        
        # get temporary index
        tmp_idx <- idx_list[[i]][r,]
        
        # get temporary variables out
        if(names(idx_list)[i] == "idx_srvjp_caa") {
          # get temporary variables out
          tmp_obs <- obs_list[[i]][,r,tmp_idx == 1,1]
          tmp_exp <- pred_list[[i]][,r,tmp_idx == 1,1]
        } else if(names(idx_list)[i] == "idx_srvdom_caa") {
          # get temporary variables out
          tmp_obs <- obs_list[[i]][,r,tmp_idx == 1,2]
          tmp_exp <- pred_list[[i]][,r,tmp_idx == 1,2]
        } else{
          tmp_obs <- obs_list[[i]][,r,tmp_idx == 1]
          tmp_exp <- pred_list[[i]][,r,tmp_idx == 1]
        }
        
        # get iss
        tmp_iss_obs <- apply(tmp_obs, 2, sum)
        tmp_obs <- t(apply(X = tmp_obs, MARGIN = 2, FUN = function(x) x / sum(x))) # normalize tmp obs 
        tmp_exp <- t(tmp_exp) # transpose temp pred variable
        
        # Set up reweighting vectors
        exp_bar <- vector() # mean expected
        obs_bar <- vector() # mean observed
        v_y <- vector() # variance
        w_denom <- vector() # weight factor in denominator
        
        for(y in 1:dim(tmp_obs)[1]) {
          exp_bar[y] <- sum(bin_list[[i]] * tmp_exp[y,]) # get mean pred comps
          obs_bar[y] <- sum(bin_list[[i]] * tmp_obs[y,]) # get mean obs comps
          v_y[y] <- sum(bin_list[[i]]^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
          w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/40) # get weights (fixing at constant 40)
        } # end y loop
        
        # Get weights 
        wts[i,r] <- 1 / var(w_denom)
        
      } # end r loop
    } # end i loop
  } # end if for spatial
  
  return(wts)
}
