#'
#' Generate tag release and recovery age-frequencies
#' due to the dimensionality this code is pretty horrendous
#' and can be quite slow, but in essence it applies the following algorithm
#' for each release event = release_year:release_region - 190 events each with a possible max_years_at_liberty x n_region recovery permuations!!! thats a 76950 different recovery age frequencies. My goodness
#' - calculate ALK
#' - calculate numbers at age of release and sex (save into a DF called release_df)
#' -- for each recovery combination - recovery year:recovery region # you can set the max number of year of recoveries you want to consider with the variable called years_at_liberty
#' -- calculate numbers at age or recovery and sex
#' --- this is done for each recovered fish!! multiply its length at release by the ALK then increment by time at liberty to get numbers at age at recovery.
#' --- link each recovery combination with a release event and save in recovery_df
# 
# source("(00) Init.R")
library(tidyverse)
library(SpatialSablefishAssessment)
library(sf)
library(here)
## survey age data
survey_age_data = readRDS(file = file.path("Data",  "Survey", "age_df_for_AF.RDS")) # created in DataCharacterisation/SurveyData.R
max_age = 31
min_age = 2
allages = min_age:max_age

length_bins = c(40, seq(from = 41.99, to = 97.99, by = 2), 200)
length_labels = paste0(c(40, seq(from = 42, to = 96, by = 2)), "-",c(seq(from = 42, to = 98, by = 2)))
length_labels = c(length_labels, "98+")
length_midpoints = seq(from=41, to=99,by=2)
## convert Age to integer
survey_age_data$age_int = as.integer(as.character(survey_age_data$Age))
## create new sex variable
survey_age_data = survey_age_data %>% mutate(Sex.x = case_when(Sex == 1  ~ "M",
                                                               Sex == 2  ~ "F",
                                                               Sex == 3  ~ "U"))

survey_age_data$length_bins = cut(survey_age_data$Length, breaks = length_bins, labels = length_midpoints)
survey_age_data$length_bins = as.numeric(as.character(survey_age_data$length_bins))
## Drop NA's because the first length bin is not a plus group then these have been assigned NA
survey_age_data = subset(survey_age_data, subset = !is.na(survey_age_data$length_bins))
## drop NAs
survey_age_data = subset(survey_age_data, subset = !is.na(survey_age_data$area_lab))
survey_age_data = survey_age_data %>% mutate(sex_length = paste0(Sex.x, length_bins))

## Bring in tag data
tag_recovery_df = readRDS(file = file.path("Data", "Tag_recoveryDF.RDS")) %>% st_drop_geometry()
tag_release_df = readRDS(file = file.path("Data", "Tag_releaseDF.RDS")) %>% st_drop_geometry()

## separate inner-gulf releases CH and CL
area_labs = data.frame(actual_lab = c("WGOA","Inner-EGOA","CGOA","WGOA", "CGOA", "EGOA","EGOA","BS","AI"),
                       obs_catch_lab = c("WGOA", "EGOA (outer & innner)","CGOA","Western Gulf of Alaska", "Central Gulf of Alaska", "West Yakutat", "East Yakutat / Southeast Alaska", "Bering Sea" , "Aleutian Islands" ))
## 
#area_labs = data.frame(actual_lab = c("WGOA","EGOA","CGOA","WGOA", "CGOA", "EGOA","EGOA","BS","AI"),
#                       obs_catch_lab = c("WGOA", "EGOA","CGOA","Western Gulf of Alaska", "Central Gulf of Alaska", "West Yakutat", "East Yakutat / Southeast Alaska", "Bering Sea" , "Aleutian Islands" ))

tag_recovery_df$area_lab = area_labs$actual_lab[match(tag_recovery_df$NAME, area_labs$obs_catch_lab)]
tag_release_df$area_lab = area_labs$actual_lab[match(tag_release_df$NAME, area_labs$obs_catch_lab)]

areas_of_interest = c("BS","AI", "WGOA", "CGOA","EGOA")
tag_recovery_df = tag_recovery_df %>% filter(area_lab %in% areas_of_interest, !is.na(area_lab))
tag_release_df = tag_release_df %>% filter(area_lab %in% areas_of_interest, !is.na(area_lab))


table(tag_recovery_df$area_lab)
table(tag_release_df$area_lab)
table(is.na(tag_recovery_df$area_lab))

## drop recoveries that we excluded in the releases
init_records = nrow(tag_recovery_df)
tag_recovery_df = tag_recovery_df %>% filter(recovery_id %in% tag_release_df$recovery_id)
cat("Recovery records dropped due to not being the release df ", init_records - nrow(tag_recovery_df), "\n")

## tag-types allowed
# - SB: SB tags were used on longline surveys from 1986-1994.
# - SA: SA tags were used on longline surveys from 1978-1987.
# - JU: released by the Japanese during the cooperative longline surveys of 1978-1983. Release and recovery data is maintained in the tag database by NMFS.
# - BK: This prefix has replaced AB and is now the main ABL tag prefix. It is used on all cruises, for juveniles as well as adults.
# - AB: This is the first tag prefix used by ABL. It was first used in 1983 and the last ones were put out in 2000. AB tags were used on juveniles as well as adults, and on inside waters cruises as well as longline surveys.
tag_types_allowed = c("AB", "BK", "JU", "SA", "SB")
table(tag_release_df$TAG_TYPE)

init_records = nrow(tag_recovery_df)
tag_recovery_df = tag_recovery_df %>% filter(TAG_TYPE.y %in% tag_types_allowed)
cat("Recovery records dropped due to tag-type ", init_records - nrow(tag_recovery_df), "\n")
init_records = nrow(tag_release_df)
tag_release_df = tag_release_df %>% filter(TAG_TYPE %in% tag_types_allowed)
cat("Release records dropped due to tag-type ", init_records - nrow(tag_recovery_df), "\n")

## subset recovery's to fixed gear i.e., longline and pt
tag_recovery_df = tag_recovery_df %>% filter(GEAR %in% c(901, 902))

## subset recoveries by the longline survey
# A recovery by the survey is identified as having an entry in the 'RCRUISE' column
# a value of 200 indicates an observer which stay as a fishery recovery
tag_recovery_df = tag_recovery_df %>% filter(is.na(RCRUISE) | RCRUISE == "200")

## prep age df
max_age = 31
min_age = 2
allages = min_age:max_age
## convert Age to integer
survey_age_data$age_int = as.integer(as.character(survey_age_data$Age))
## create new sex variable

survey_age_data = survey_age_data %>% mutate(Sex.x = case_when(Sex == 1  ~ "M",
                                                               Sex == 2  ~ "F",
                                                               Sex == 3  ~ "U"))

length_bins = c(40, seq(from = 41.99, to = 97.99, by = 2), 200)
length_labels = paste0(c(40, seq(from = 42, to = 96, by = 2)), "-",c(seq(from = 42, to = 98, by = 2)))
length_labels = c(length_labels, "98+")
length_midpoints = seq(from=41, to=99,by=2)

survey_age_data$length_bins = cut(survey_age_data$Length, breaks = length_bins, labels = length_midpoints, ordered_result  = T)
tag_recovery_df$recovery_length_bin = cut(tag_recovery_df$RSIZE/10, breaks = length_bins, labels = length_midpoints, ordered_result  = T)
tag_recovery_df$release_length_bins = cut(tag_recovery_df$HSIZE.x/10, breaks = length_bins, labels = length_midpoints, ordered_result  = T)
tag_release_df$release_length_bins = cut(tag_release_df$HSIZE/10, breaks = length_bins, labels = length_midpoints, ordered_result  = T)
## years at liberty
tag_recovery_df$year_at_liberty = as.numeric(tag_recovery_df$recovery_year) - as.numeric(tag_recovery_df$release_year)

sex_age_lvls = paste0(rep(c("M","F"), each = length(allages)), allages)
age_for_report = rep(allages, 2)
sex_for_report= rep(c("M","F"), each = length(allages))
survey_age_data = survey_age_data %>% mutate(sex_age = paste0(Sex.x, age_int))
## drop entries that are not in our sex combos
survey_age_data = survey_age_data %>% filter(sex_age %in% sex_age_lvls)
survey_age_data$sex_age = factor(survey_age_data$sex_age, levels = sex_age_lvls, ordered = T)

## years with age-date
years_with_age_data = sort(unique(survey_age_data$Year.x))

############################
## Visualise 
## empirical recovery movement matrix
##############################
temp = tag_recovery_df %>% inner_join(tag_release_df, by = "recovery_id")
## when we merge area_lab.x is recovery region and area_lab.y is release region
temp_release = tag_release_df %>% full_join(tag_recovery_df, by = "recovery_id")
## when we merge area_lab.x is RELEASE region and area_lab.y is RECOVERY region
temp_release = temp_release %>% group_by(area_lab.y, area_lab.x) %>% summarise(releases = n(), recoveries = sum(!is.na(TAG_TYPE.x)))
temp_release$area_lab.y = factor(temp_release$area_lab.y, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
temp_release$area_lab.x = factor(temp_release$area_lab.x, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
sum_tag_release_df = temp_release %>% group_by(area_lab.x) %>% mutate(prop = recoveries / sum(recoveries), releases= sum(releases))
# drop NA regions
sum_tag_release_df = sum_tag_release_df %>% filter(!is.na(area_lab.y))

sum_tag_release_df$area_lab_release = paste0(sum_tag_release_df$area_lab.x, "\n(", sum_tag_release_df$releases,")")
unique(sum_tag_release_df$area_lab_release)
sum_tag_release_df$area_lab_release = factor(sum_tag_release_df$area_lab_release, levels = rev(c(    "BS\n(28126)" , "AI\n(18917)","WGOA\n(30875)",   "CGOA\n(81164)",  "EGOA\n(117342)" )), ordered = T)

ggplot(sum_tag_release_df, aes(x = area_lab.y, y = area_lab_release, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(x = area_lab.y, y = area_lab_release, label = round(recoveries,2)), color = "black", size = 4) +
  labs(x = "Recovery region", y = "Release region (Number of releases)", fill = "Proportion of\nrecoveries") +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14)
  )
ggsave(filename = file.path(DIR$input_figs_5A, "empirical_tag_recoveries_with_releases.png"), width = 7, height = 6)

temp_recovery = temp %>% group_by(area_lab.y, area_lab.x) %>% summarise(recoveries = n())
temp_recovery$area_lab.y = factor(temp_recovery$area_lab.y, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
temp_recovery$area_lab.x = factor(temp_recovery$area_lab.x, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
sum_tag_recap_df = temp_recovery %>% group_by(area_lab.y) %>% mutate(prop = recoveries / sum(recoveries))
ggplot(sum_tag_recap_df, aes(x = area_lab.x, y = area_lab.y, fill = recoveries)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(x = area_lab.x, y = area_lab.y, label = round(prop,2)), color = "black", size = 4) +
  labs(x = "Recovery region", y = "Release region", fill = "Number of\nrecoveries")
ggsave(filename = file.path(DIR$input_figs_5A, "empirical_tag_recoveries.png"), width = 7, height = 6)



## redo the recapture probabilities by Dana's length groups
table(is.na(temp$length_group))
# drop NA length groups
temp = temp %>% filter(!is.na(length_group))
# Recode factor levels by name
temp$length_group = factor(temp$length_group)
temp$length_group_f = recode_factor(temp$length_group, small  = "Small (<56.5 cm)", medium  = "Medium (56.5-66.5 cm)", large  = "Large (>66.5 cm)")
temp_recovery = temp %>% group_by(area_lab.y, area_lab.x, length_group_f) %>% summarise(recoveries = n())
temp_recovery$area_lab.y = factor(temp_recovery$area_lab.y, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
temp_recovery$area_lab.x = factor(temp_recovery$area_lab.x, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
sum_tag_recap_df = temp_recovery %>% group_by(area_lab.y,length_group_f) %>% mutate(prop = recoveries / sum(recoveries))
ggplot(sum_tag_recap_df, aes(x = area_lab.x, y = area_lab.y, fill = recoveries)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(x = area_lab.x, y = area_lab.y, label = round(prop,2)), color = "black", size = 4) +
  labs(x = "Recovery region", y = "Release region", fill = "Recoveries") +
  facet_wrap(~length_group_f)
ggsave(filename = file.path(DIR$input_figs_5A, "empirical_tag_recoveries_by_length_group.png"), width = 12, height = 6)

## when we merge area_lab.x is recovery region and area_lab.y is release region
temp_release = tag_release_df %>% full_join(tag_recovery_df, by = "recovery_id")

temp_release$length_group = factor(temp_release$length_group)
temp_release$length_group_f = recode_factor(temp_release$length_group, small  = "Small (<56.5 cm)", medium  = "Medium (56.5-66.5 cm)", large  = "Large (>66.5 cm)")
## when we merge area_lab.x is RELEASE region and area_lab.y is RECOVERY region
temp_release = temp_release %>% group_by(area_lab.y, area_lab.x, length_group_f) %>% summarise(releases = n(), recoveries = sum(!is.na(TAG_TYPE.x)))
temp_release$area_lab.y = factor(temp_release$area_lab.y, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
temp_release$area_lab.x = factor(temp_release$area_lab.x, levels = rev(c("BS", "AI","WGOA","CGOA","EGOA")), ordered = T)
sum_tag_release_df = temp_release %>% group_by(area_lab.x, length_group_f) %>% mutate(prop = recoveries / sum(recoveries), releases= sum(releases))
# drop NA regions
sum_tag_release_df = sum_tag_release_df %>% filter(!is.na(area_lab.y))
# drop NA length groups
sum_tag_release_df = sum_tag_release_df %>% filter(!is.na(length_group_f))

sum_tag_release_df$area_lab_release = paste0(sum_tag_release_df$area_lab.x, "\n(", sum_tag_release_df$releases,")")
unique(sum_tag_release_df$area_lab_release)
#sum_tag_release_df$area_lab_release = factor(sum_tag_release_df$area_lab_release, levels = rev(c(    "BS\n(28126)" , "AI\n(18917)","WGOA\n(30875)",   "CGOA\n(81164)",  "EGOA\n(117342)" )), ordered = T)

ggplot(sum_tag_release_df, aes(x = area_lab.y, y = area_lab.x, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(x = area_lab.y, y = area_lab.x, label = round(recoveries,2)), color = "black", size = 4) +
  labs(x = "Recovery region", y = "Release region (Number of releases)", fill = "Proportion of\nrecoveries") +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14)
  ) +
  facet_wrap(~length_group_f)

ggsave(filename = file.path(DIR$input_figs_5A, "empirical_tag_recoveries_by_length_group.png"), width = 12, height = 6)

## years with release
years_with_release_events = as.numeric(sort(unique(tag_release_df$release_year)))
n_release_years = length(years_with_release_events)
regions = unique(tag_release_df$area_lab)
years_at_liberty = 15;

## age_fish - used to calculate age of at recovery
age_this_fish = function(nums_at_age, increment_years, allages = 2:31) {
  copy_nage = nums_at_age
  tmp = rep(0, length(nums_at_age))
  tmp[(increment_years+1):length(allages)] = copy_nage[1:(length(allages) - increment_years)]
  tmp[length(allages)] = sum(copy_nage[(length(allages) - increment_years):length(allages)])
  
  tmp[length(allages) + (increment_years+1):length(allages)] = copy_nage[length(allages) + 1:(length(allages) - increment_years)]
  tmp[length(allages) + length(allages)] = sum(copy_nage[length(allages) +(length(allages) - increment_years):length(allages)])
  return(tmp)
}

## create a global ALK used to impute length bins that 
## 
global_freq_len_age = xtabs(~sex_length + sex_age, data = survey_age_data, drop.unused.levels = FALSE, addNA = TRUE)
## are there any length bins with no data
global_rowSums = rowSums(global_freq_len_age)
if(any(global_rowSums == 0)){
  global_rowSums[global_rowSums == 0]  = 1 ## this is to remove NaNs
}
global_ALK_year = sweep(global_freq_len_age, 1,STATS = global_rowSums, FUN = "/")


## calculate Age-frequency for each region and year
## simialr to ALK AF method
## calculate ALK for a year across all regions.
tag_release_event = 1;
release_df = NULL
recovery_df = NULL
temp_age_vector = NULL
for(y_ndx in 1:n_release_years) {
  if(years_with_release_events[y_ndx] < years_with_age_data[1]) {
    this_pooled_age_df = survey_age_data %>% filter(Year.x %in% years_with_age_data[1:2])
  } else {
    if(years_with_release_events[y_ndx] %in% years_with_age_data) {
      this_pooled_age_df = survey_age_data %>% filter(Year.x == years_with_release_events[y_ndx])
    } else {
      ## grab the age-length key in the closest year 
      this_pooled_age_df = survey_age_data %>% filter(Year.x == years_with_age_data[which.min(abs(years_with_release_events[y_ndx] - years_with_age_data))])
    }
  }
  ## build age length key
  pooled_freq_len_age = xtabs(~length_bins + sex_age, data = this_pooled_age_df, drop.unused.levels = FALSE, addNA = F)
  pooled_row_totals = rowSums(pooled_freq_len_age)
  # What to do with lengths with zeros
  zero_ndx = which(pooled_row_totals == 0)
  pooled_row_totals[pooled_row_totals == 0] = 1
  pooled_ALK_year = sweep(pooled_freq_len_age, 1,STATS = pooled_row_totals, FUN = "/")
  pooled_ALK_year[zero_ndx, ] = global_ALK_year[zero_ndx, ]
  
  pooled_num_at_len_A = this_pooled_age_df %>%
    group_by(length_bins) %>%
    summarise(NUMBERS = n()) %>%
    complete(length_bins, fill = list(NUMBERS = 0))
  
  for(region_ndx in 1:length(regions)) {
    releases_this_year_region = tag_release_df %>% filter(release_year == years_with_release_events[y_ndx], area_lab == regions[region_ndx], !is.na(release_length_bins))
    
    if(nrow(releases_this_year_region) > 0) {
      
      num_at_len = releases_this_year_region %>%
        group_by(release_length_bins) %>%
        summarise(NUMBERS = n()) %>%
        complete(release_length_bins, fill = list(NUMBERS = 0))
      ## match scaled Length bins with age-length key
      num_at_len = num_at_len %>% filter(release_length_bins %in% pooled_num_at_len_A$length_bins)
      ## get ALK age-dist
      N_a_release = apply(pooled_ALK_year, 2, function(x) {sum(x * num_at_len$NUMBERS, na.rm = T)})
      
      ## save in a release DF
      tmp_release_df = data.frame(release_event_id = tag_release_event, Nage_at_release = N_a_release, release_year = years_with_release_events[y_ndx], region_release = regions[region_ndx],
                                  age = age_for_report, sex = sex_for_report, sex_age = sex_age_lvls)
      release_df = rbind(release_df, tmp_release_df)
      
      ## check if we lost some fish
      if(round(sum(N_a_release)) != round(sum(num_at_len$NUMBERS))) {
        print(paste0("we lost ", sum(num_at_len$NUMBERS) - sum(N_a_release), " fish in year ", years_with_release_events[y_ndx], " and region = ", regions[region_ndx], " during ALK"))
      }
      ## get the recoveries from this release event up to 
      related_recoveries = tag_recovery_df %>% filter(recovery_id %in% releases_this_year_region$recovery_id)
      ## only want recoveries caught within a certain amount of years at liberty
      related_recoveries = related_recoveries %>% filter(year_at_liberty <= years_at_liberty)
      ## cut these up by region and calculate the corresponding AF
      if(nrow(related_recoveries) > 0) {
        print(paste0("In year ", years_with_release_events[y_ndx], " and region = ", regions[region_ndx], " found ", nrow(related_recoveries) , " recoveries of ", nrow(releases_this_year_region), " releases"))
        ## iterate over year then region
        for(recovery_year_ndx in 0:years_at_liberty) {
          for(recovery_region_ndx in 1:length(regions)) {
            region_recoveries = related_recoveries %>% filter(area_lab  == regions[recovery_region_ndx], year_at_liberty == recovery_year_ndx)
            if(nrow(region_recoveries) > 0) {
              N_a_recovery = rep(0, length(allages) * 2)
              for(recovered_fish in 1:nrow(region_recoveries)) {
                this_fish = region_recoveries[recovered_fish, ]
                length_dist = this_fish %>% mutate(length_ndx = n()) %>% dplyr::select(length_ndx, release_length_bins) %>% complete(release_length_bins, fill = list(length_ndx = 0))
                ## this will be a distribution
                age_and_sex_at_release = apply(pooled_ALK_year, 2, function(x) {sum(x * length_dist$length_ndx, na.rm = T)})
                ## increment age
                incr_age = as.numeric(this_fish$recovery_year) - as.numeric(this_fish$release_year)
                age_and_sex_at_recovery = age_this_fish(age_and_sex_at_release,incr_age) 
                N_a_recovery = N_a_recovery + age_and_sex_at_recovery 
              }
              ## save in a release DF
              tmp_recovery_df = data.frame(release_event_id = tag_release_event, Nage_at_recovery = N_a_recovery, release_year = years_with_release_events[y_ndx], region_release = regions[region_ndx], recovery_region = regions[recovery_region_ndx], recovery_year = as.numeric(this_fish$recovery_year),
                                           age = age_for_report, sex = sex_for_report, sex_age = sex_age_lvls)
              recovery_df = rbind(recovery_df, tmp_recovery_df)
            }
          }
        }
      }
      tag_release_event = tag_release_event + 1;
    }
  }
}

cat("total tags released = ", sum(recovery_df$Nage_at_recovery), "\n")
cat("total tags recovered = ", sum(release_df$Nage_at_release), "\n")


saveRDS(recovery_df, file = file.path(DIR$data5A, "Tag_recovery_summarised.RDS"))
saveRDS(release_df, file = file.path(DIR$data5A, "Tag_release_summarised.RDS"))

## do some simple arithmetic


