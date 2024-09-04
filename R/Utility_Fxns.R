

#' Title
#'
#' @param sablefish_control_filename Input file path for sablefish .ctl file
#'
#' @return List of ctl variables
#' @export
#'
#' @examples
digest_sab_ctl = function(sablefish_control_filename) {
  ctl_files = readLines(sablefish_control_filename) 
  ## only get the key things
  ## rec start year
  rec_styr <-as.numeric(unlist(strsplit(ctl_files[grep("styr for rec devs est",ctl_files)+1],split=" "))) #Get mean recruitment used for B40 calc from report file
  ## end yr for recrum
  rec_endyr <- as.numeric(unlist(strsplit(ctl_files[grep("end yr for rec devs est",ctl_files)+1],split=" "))) #Get mean recruitment used for B40 calc from report file
  ## SR type = 
  SRtype = as.numeric(unlist(strsplit(ctl_files[grep("#SR type",ctl_files)+1],split=" "))) #Get mean recruitment used for B40 calc from report file
  
  return(list(rec_styr = rec_styr, rec_endyr = rec_endyr, SRtype = SRtype)) 
}


#' Title
#'
#' @param sablefish_input_filename Input file path for sablefish data inputs
#'
#' @return List of sablefish .dat file inputs
#' @export
#'
#' @examples
digest_sab_input = function(sablefish_input_filename) {
  input_files = readLines(sablefish_input_filename) 
  yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("Start and end years",input_files)+1],split=" ")))[1:2])
  recage = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("Start and end years",input_files)+2],split=" ")))[1])
  nlenbins = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("Start and end years",input_files)+3],split=" ")))[1]) 
  length_bins = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# length bins",input_files)+1],split=" "))))
  spawn_month = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# length bins",input_files)+2],split=" ")))[1])
  male_60_96_ndx = grep(pattern ="# Size-at-age, 1960-1996",input_files)[1]
  female_60_96_ndx = grep(pattern ="# Size-at-age, 1960-1996",input_files)[2]
  male_length_at_age_60_96 = matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(male_60_96_ndx + 1):(male_60_96_ndx + length(length_bins))],split=" ")))), byrow = T, ncol = length(length_bins))
  female_length_at_age_60_96 = matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(female_60_96_ndx + 1):(female_60_96_ndx + length(length_bins))],split=" ")))), byrow = T, ncol = length(length_bins))
  male_97_22_ndx = grep(pattern ="# Size-at-age, 1997-2021",input_files)[1]
  female_97_22_ndx = grep(pattern ="# Size-at-age, 1997-2021",input_files)[2]
  male_length_at_age_97_22 = matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(male_97_22_ndx + 1):(male_97_22_ndx + length(length_bins))],split=" ")))), byrow = T, ncol = length(length_bins))
  female_length_at_age_97_22 = matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(female_97_22_ndx + 1):(female_97_22_ndx + length(length_bins))],split=" ")))), byrow = T, ncol = length(length_bins))
  
  jap_fishery_ll_age_length_ndx = grep(pattern ="# Size-at-age, historical, JPN Length Data, Unsexed",input_files)[1]
  jap_fishery_ll_age_length_transition =matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(jap_fishery_ll_age_length_ndx + 1):(jap_fishery_ll_age_length_ndx + length(length_bins))],split=" ")))), byrow = T, ncol = length(length_bins)) 
  
  age_error_ndx = grep(pattern ="#Age age transition matrix",input_files)
  ageing_error = matrix(suppressWarnings(as.numeric(unlist(strsplit(input_files[(age_error_ndx + 1):(age_error_ndx + nrow(female_length_at_age_97_22))],split=" ")))), byrow = T, nrow = nrow(female_length_at_age_97_22))
  
  srv_dom_ll_bio_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Domestic longline survey RPN",input_files)+5],split=" "))))
  srv_dom_ll_bio_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Domestic longline survey RPN",input_files)+7)],split=" "))))
  srv_dom_ll_bio_se = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Domestic longline survey RPN",input_files)+9)],split=" "))))
  
  srv_jap_ll_bio_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese/coop longline survey RPN",input_files)+5],split=" "))))
  srv_jap_ll_bio_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Japanese/coop longline survey RPN",input_files)+7)],split=" "))))
  srv_jap_ll_bio_se = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Japanese/coop longline survey RPN",input_files)+9)],split=" "))))
  
  ll_cpue_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Domestic Longline Fishery CPUE",input_files)+5],split=" "))))
  ll_cpue_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Domestic Longline Fishery CPUE",input_files)+7)],split=" "))))
  ll_cpue_se = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Domestic Longline Fishery CPUE",input_files)+9)],split=" "))))
  
  ll_age_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Fixed Gear Fishery Age Composition",input_files)+5],split=" "))))
  ll_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Fixed Gear Fishery Age Composition",input_files)+7],split=" "))))
  ll_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Fixed Gear Fishery Age Composition",input_files)+9):(grep("# Fixed Gear Fishery Age Composition",input_files) +9 + length(ll_age_yrs))],split=" "))))
  ll_obs = subset(ll_obs, subset = !is.na(ll_obs))
  ll_obs_mat = matrix(ll_obs, nrow = length(ll_age_yrs), byrow = T)
  
  trwl_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Trawl gear fishery length compositions",input_files)+5],split=" "))))
  trwl_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Trawl gear fishery length compositions",input_files)+7],split=" "))))
  trwl_lgth_m_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("Observed trawl fishery size compositions MALE",input_files)+1):(grep("Observed trawl fishery size compositions MALE",input_files) + length(trwl_lgth_yrs))],split=" "))))
  trwl_lgth_m_obs = subset(trwl_lgth_m_obs, subset = !is.na(trwl_lgth_m_obs))
  trwl_lgth_m_obs_mat = matrix(trwl_lgth_m_obs, nrow = length(trwl_lgth_yrs), byrow = T)
  trwl_lgth_f_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed trawl fishery size compositions Females",input_files)+1):(grep("# Observed trawl fishery size compositions Females",input_files) + length(trwl_lgth_yrs))],split=" "))))
  trwl_lgth_f_obs = subset(trwl_lgth_f_obs, subset = !is.na(trwl_lgth_f_obs))
  trwl_lgth_f_obs_mat = matrix(trwl_lgth_f_obs, nrow = length(trwl_lgth_yrs), byrow = T)
  
  
  srv_dom_ll_age_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Domestic LL survey Age Composition",input_files)+5],split=" "))))
  srv_dom_ll_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Domestic LL survey Age Composition",input_files)+7],split=" "))))
  srv_dom_ll_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Domestic LL survey Age Composition",input_files)+9):(grep("# Domestic LL survey Age Composition",input_files) +9 + length(srv_dom_ll_age_yrs))],split=" "))))
  srv_dom_ll_obs = subset(srv_dom_ll_obs, subset = !is.na(srv_dom_ll_obs))
  srv_dom_ll_obs_mat = matrix(srv_dom_ll_obs, nrow = length(srv_dom_ll_age_yrs), byrow = T)
  
  srv_dom_ll_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Domestic LL Survey Size Composition",input_files)+5],split=" "))))
  srv_dom_ll_lgth_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Domestic LL Survey Size Composition",input_files)+7],split=" "))))
  srv_dom_ll_lgth_m_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed LL survey size compositions MALE",input_files)+1):(grep("# Observed LL survey size compositions MALE",input_files) + length(srv_dom_ll_lgth_yrs))],split=" "))))
  srv_dom_ll_lgth_m_obs = subset(srv_dom_ll_lgth_m_obs, subset = !is.na(srv_dom_ll_lgth_m_obs))
  srv_dom_ll_lgth_m_obs_mat = matrix(srv_dom_ll_lgth_m_obs, nrow = length(srv_dom_ll_lgth_yrs), byrow = T)
  srv_dom_ll_lgth_f_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed LL survey size compositions FEMALE",input_files)+1):(grep("# Observed LL survey size compositions FEMALE",input_files) + length(srv_dom_ll_lgth_yrs))],split=" "))))
  srv_dom_ll_lgth_f_obs = subset(srv_dom_ll_lgth_f_obs, subset = !is.na(srv_dom_ll_lgth_f_obs))
  srv_dom_ll_lgth_f_obs_mat = matrix(srv_dom_ll_lgth_f_obs, nrow = length(srv_dom_ll_lgth_yrs), byrow = T)
  
  srv_jap_ll_age_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese LL survey Age Composition",input_files)+5],split=" "))))
  srv_jap_ll_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese LL survey Age Composition",input_files)+7],split=" "))))
  srv_jap_ll_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Japanese LL survey Age Composition",input_files)+9):(grep("# Japanese LL survey Age Composition",input_files) +9 + length(srv_jap_ll_age_yrs))],split=" "))))
  srv_jap_ll_obs = subset(srv_jap_ll_obs, subset = !is.na(srv_jap_ll_obs))
  srv_jap_ll_obs_mat = matrix(srv_jap_ll_obs, nrow = length(srv_jap_ll_age_yrs), byrow = T)
  
  srv_jap_ll_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese LL survey size Composition",input_files)+5],split=" "))))
  srv_jap_ll_lgth_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese LL survey size Composition",input_files)+7],split=" "))))
  srv_jap_ll_lgth_m_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed JPN LL survey size compositions Male",input_files)+1):(grep("# Observed JPN LL survey size compositions Male",input_files) + length(srv_jap_ll_lgth_yrs))],split=" "))))
  srv_jap_ll_lgth_m_obs = subset(srv_jap_ll_lgth_m_obs, subset = !is.na(srv_jap_ll_lgth_m_obs))
  srv_jap_ll_lgth_m_obs_mat = matrix(srv_jap_ll_lgth_m_obs, nrow = length(srv_jap_ll_lgth_yrs), byrow = T)
  srv_jap_ll_lgth_f_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed JPN LL survey size compositions FeMale",input_files)+1):(grep("# Observed JPN LL survey size compositions FeMale",input_files) + length(srv_jap_ll_lgth_yrs))],split=" "))))
  srv_jap_ll_lgth_f_obs = subset(srv_jap_ll_lgth_f_obs, subset = !is.na(srv_jap_ll_lgth_f_obs))
  srv_jap_ll_lgth_f_obs_mat = matrix(srv_jap_ll_lgth_f_obs, nrow = length(srv_jap_ll_lgth_yrs), byrow = T)
  
  srv_nmfs_trwl_age_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# GOA trawl survey historical age composition",input_files)+5],split=" "))))
  srv_nmfs_trwl_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# GOA trawl survey historical age composition",input_files)+7],split=" "))))
  srv_nmfs_trwl_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# GOA trawl survey historical age composition",input_files)+9):(grep("# GOA trawl survey historical age composition",input_files) +9 + length(srv_nmfs_trwl_age_yrs))],split=" "))))
  srv_nmfs_trwl_obs = subset(srv_nmfs_trwl_obs, subset = !is.na(srv_nmfs_trwl_obs))
  srv_nmfs_trwl_obs_mat = matrix(srv_nmfs_trwl_obs, nrow = length(srv_nmfs_trwl_age_yrs), byrow = T)
  
  srv_nmfs_trwl_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# GOA trawl survey size composition",input_files)+5],split=" "))))
  srv_nmfs_trwl_lgth_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# GOA trawl survey size composition",input_files)+7],split=" "))))
  srv_nmfs_trwl_lgth_m_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed Trawl Survey size comps Male",input_files)+1):(grep("# Observed Trawl Survey size comps Male",input_files) + length(srv_nmfs_trwl_lgth_yrs))],split=" "))))
  srv_nmfs_trwl_lgth_m_obs = subset(srv_nmfs_trwl_lgth_m_obs, subset = !is.na(srv_nmfs_trwl_lgth_m_obs))
  srv_nmfs_trwl_lgth_m_obs_mat = matrix(srv_nmfs_trwl_lgth_m_obs, nrow = length(srv_nmfs_trwl_lgth_yrs), byrow = T)
  srv_nmfs_trwl_lgth_f_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed Trawl Survey size comps Females",input_files)+1):(grep("# Observed Trawl Survey size comps Females",input_files) + length(srv_nmfs_trwl_lgth_yrs))],split=" "))))
  srv_nmfs_trwl_lgth_f_obs = subset(srv_nmfs_trwl_lgth_f_obs, subset = !is.na(srv_nmfs_trwl_lgth_f_obs))
  srv_nmfs_trwl_lgth_f_obs_mat = matrix(srv_nmfs_trwl_lgth_f_obs, nrow = length(srv_nmfs_trwl_lgth_yrs), byrow = T)
  
  ll_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Fixed gear fishery length compositions",input_files)+5],split=" "))))
  ll_lgth_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# U.S. Fixed gear fishery length compositions",input_files)+7],split=" "))))
  ll_lgth_m_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed LL fishery size compositions MALE",input_files)+1):(grep("# Observed LL fishery size compositions MALE",input_files) + length(ll_lgth_yrs))],split=" "))))
  ll_lgth_m_obs = subset(ll_lgth_m_obs, subset = !is.na(ll_lgth_m_obs))
  ll_lgth_m_obs_mat = matrix(ll_lgth_m_obs, nrow = length(ll_lgth_yrs), byrow = T)
  ll_lgth_f_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed LL fishery size compositions Females",input_files)+1):(grep("# Observed LL fishery size compositions Females",input_files) + length(ll_lgth_yrs))],split=" "))))
  ll_lgth_f_obs = subset(ll_lgth_f_obs, subset = !is.na(ll_lgth_f_obs))
  ll_lgth_f_obs_mat = matrix(ll_lgth_f_obs, nrow = length(ll_lgth_yrs), byrow = T)
  
  srv_nmfs_trwl_bio_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# GOA Trawl Survey Biomass",input_files)+5],split=" "))))
  srv_nmfs_trwl_bio_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# GOA Trawl Survey Biomass",input_files)+7)],split=" "))))
  srv_nmfs_trwl_bio_se = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# GOA Trawl Survey Biomass",input_files)+9)],split=" "))))
  
  srv_jap_fishery_ll_bio_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese longline fishery CPUE RPW",input_files)+5],split=" "))))
  srv_jap_fishery_ll_bio_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Japanese longline fishery CPUE RPW",input_files)+7)],split=" "))))
  srv_jap_fishery_ll_bio_se = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Japanese longline fishery CPUE RPW",input_files)+9)],split=" "))))
  
  srv_jap_fishery_ll_lgth_yrs = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese longline fishery unsexed lengths",input_files)+5],split=" "))))
  srv_jap_fishery_ll_lgth_sample_size = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Japanese longline fishery unsexed lengths",input_files)+7],split=" "))))
  srv_jap_fishery_ll_lgth_obs = suppressWarnings(as.numeric(unlist(strsplit(input_files[(grep("# Observed JPN LL fishery size compositions",input_files)+1):(grep("# Observed JPN LL fishery size compositions",input_files) + length(srv_jap_fishery_ll_lgth_yrs))],split=" "))))
  srv_jap_fishery_ll_lgth_obs = subset(srv_jap_fishery_ll_lgth_obs, subset = !is.na(srv_jap_fishery_ll_lgth_obs))
  srv_jap_fishery_ll_lgth_obs_mat = matrix(srv_jap_fishery_ll_lgth_obs, nrow = length(srv_jap_fishery_ll_lgth_yrs), byrow = T)
  
  ll_Catch = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Fixed gear fishery catch",input_files)+1],split=" "))))
  trwl_Catch = suppressWarnings(as.numeric(unlist(strsplit(input_files[grep("# Trawl gear fishery catch",input_files)+1],split=" "))))
  
  
  
  return(list(spawn_month = spawn_month, length_bins = length_bins, nlenbins = nlenbins, recage = recage, yrs = yrs, male_length_at_age_60_96= male_length_at_age_60_96, 
              female_length_at_age_60_96 = female_length_at_age_60_96, male_length_at_age_97_22 = male_length_at_age_97_22, female_length_at_age_97_22 = female_length_at_age_97_22, ageing_error = ageing_error,
              ll_age_yrs = ll_age_yrs, ll_sample_size = ll_sample_size, ll_obs_mat = ll_obs_mat, jap_fishery_ll_age_length_transition = jap_fishery_ll_age_length_transition,
              trwl_lgth_yrs = trwl_lgth_yrs, trwl_sample_size = trwl_sample_size, trwl_lgth_f_obs_mat = trwl_lgth_f_obs_mat, trwl_lgth_m_obs_mat = trwl_lgth_m_obs_mat,
              srv_dom_ll_bio_yrs = srv_dom_ll_bio_yrs, srv_dom_ll_obs = srv_dom_ll_bio_obs, srv_dom_ll_se = srv_dom_ll_bio_se,
              srv_jap_ll_bio_yrs = srv_jap_ll_bio_yrs, srv_jap_ll_bio_obs = srv_jap_ll_bio_obs, srv_jap_ll_bio_se = srv_jap_ll_bio_se,
              ll_cpue_yrs = ll_cpue_yrs, ll_cpue_obs = ll_cpue_obs, ll_cpue_se = ll_cpue_se,
              srv_dom_ll_age_yrs = srv_dom_ll_age_yrs, srv_dom_ll_sample_size = srv_dom_ll_sample_size, srv_dom_ll_obs_mat = srv_dom_ll_obs_mat,
              srv_dom_ll_lgth_yrs = srv_dom_ll_lgth_yrs, srv_dom_ll_lgth_sample_size = srv_dom_ll_lgth_sample_size, srv_dom_ll_lgth_m_obs_mat = srv_dom_ll_lgth_m_obs_mat, srv_dom_ll_lgth_f_obs_mat= srv_dom_ll_lgth_f_obs_mat,
              srv_jap_ll_age_yrs = srv_jap_ll_age_yrs, srv_jap_ll_sample_size = srv_jap_ll_sample_size, srv_jap_ll_obs_mat = srv_jap_ll_obs_mat,
              srv_jap_ll_lgth_yrs = srv_jap_ll_lgth_yrs, srv_jap_ll_lgth_sample_size = srv_jap_ll_lgth_sample_size, srv_jap_ll_lgth_m_obs_mat = srv_jap_ll_lgth_m_obs_mat, srv_jap_ll_lgth_f_obs_mat= srv_jap_ll_lgth_f_obs_mat,
              srv_nmfs_trwl_age_yrs = srv_nmfs_trwl_age_yrs, srv_nmfs_trwl_sample_size = srv_nmfs_trwl_sample_size, srv_nmfs_trwl_obs_mat = srv_nmfs_trwl_obs_mat,
              srv_nmfs_trwl_lgth_yrs = srv_nmfs_trwl_lgth_yrs, srv_nmfs_trwl_lgth_sample_size = srv_nmfs_trwl_lgth_sample_size, srv_nmfs_trwl_lgth_m_obs_mat = srv_nmfs_trwl_lgth_m_obs_mat, srv_nmfs_trwl_lgth_f_obs_mat= srv_nmfs_trwl_lgth_f_obs_mat,
              ll_lgth_yrs = ll_lgth_yrs, ll_lgth_sample_size = ll_lgth_sample_size, ll_lgth_m_obs_mat = ll_lgth_m_obs_mat, ll_lgth_f_obs_mat= ll_lgth_f_obs_mat,
              srv_nmfs_trwl_bio_yrs = srv_nmfs_trwl_bio_yrs, srv_nmfs_trwl_obs = srv_nmfs_trwl_bio_obs, srv_nmfs_trwl_se = srv_nmfs_trwl_bio_se,
              srv_jap_fishery_ll_bio_yrs = srv_jap_fishery_ll_bio_yrs, srv_jap_fishery_ll_obs = srv_jap_fishery_ll_bio_obs, srv_jap_fishery_ll_se = srv_jap_fishery_ll_bio_se,
              srv_jap_fishery_ll_lgth_yrs = srv_jap_fishery_ll_lgth_yrs, srv_jap_fishery_ll_lgth_sample_size = srv_jap_fishery_ll_lgth_sample_size, srv_jap_fishery_ll_lgth_obs_mat = srv_jap_fishery_ll_lgth_obs_mat,
              ll_Catch = ll_Catch, trwl_Catch = trwl_Catch
              
  )) 
  
}

# Indicator function
indicator_fun = function(x) {
  ifelse(is.na(sum(x)), 0, 1)
}