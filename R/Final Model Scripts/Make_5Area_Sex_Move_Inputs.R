# Purpose: To make data inputs for 5-area age based and sex-based movement models 
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 9/12/24

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

scenario = "1960_start" # Model scenario
fig_path = here("Figs", "5-Area-1960")
out_path = here("Output", "Final Models", "5-Area-1960")

# Read in tag data
tag_recovery_df = readRDS(file = here("Data", "5-Area", "Tag_recovery_summarised.RDS"))
tag_release_df = readRDS(file = here("Data", "5-Area", "Tag_release_summarised.RDS"))

# Load in data and parameters from 5 area base models wihtout age-based movement
data = readRDS(file = here(out_path, "data.RDS"))
parameters = readRDS(file = here(out_path, "parameters.RDS"))
multiple_shoot = T # whether to do multiple shooting
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))

# Other specifications
n_sexes = 2 # number of sexes
data$age_based_movement = 0 # age based movement

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
data$obs_tag_recovery = array(0, dim = c(n_sexes, data$n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), data$n_regions, length(tag_recovery_years)))
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
            data$obs_tag_recovery[,model_tag_release_ndx, recovery_region_ndx, y_ndx] = ((specific_recovery_df %>% group_by(sex) %>% summarise(Nage_at_recovery = sum(Nage_at_recovery)))$Nage_at_recovery)
          }
        }
      }
    }
  }
}

# Set up movement parameterization ----------------------------------------

data$n_movement_age_blocks = 1 # three age blocks
data$movement_age_block_indicator = c(rep(0, 5), rep(0, length(6:12)), rep(0, length(13:30))) 
data$n_movement_sex_blocks = 2 # sex blocks
data$movement_sex_block_indicator = c(0,1)

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
move_matrix = array(0, dim = c(data$n_regions, data$n_regions)) # temporary matrix for feeding initial values
data$fixed_movement_matrix = array(0, dim = c(data$n_regions, data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks, data$n_movement_sex_blocks))
# set diagnoals equal to 1
# diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1]) = diag(data$fixed_movement_matrix[,,1,2]) = 1
diag(move_matrix) = diag(data$fixed_movement_matrix[,,1,1,1])  = 1
move_matrix = move_matrix + rlnorm(n = data$n_regions * data$n_regions, log(0.01), 0.5) # add random normal draws to make no equal to 1
# renormalise
move_matrix = sweep(move_matrix, 1, STATS = rowSums(move_matrix), "/")

# set up movement pars
parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks, data$n_movement_age_blocks, data$n_movement_sex_blocks))
# parameters$transformed_movement_pars = array(NA, c(data$n_regions - 1, ncol = data$n_regions, data$n_movement_time_blocks))

for(i in 1:data$n_regions) {
  for(t in 1:data$n_movement_time_blocks) {
    for(a in 1:data$n_movement_age_blocks) {
      for(s in 1:data$n_movement_sex_blocks) {
        parameters$transformed_movement_pars[,i,t,a,s] = simplex(move_matrix[i,])
      }
    } # end a block
  } # end t block
} # end i region

# Report out data and parameters 
saveRDS(data, file.path(out_path, "data_5area_agesexmove.RDS"))
saveRDS(parameters, file.path(out_path, "parameters_5area_agesexmove.RDS"))

