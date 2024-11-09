# Purpose: To plot high reoslution data explorations conducted for manuscript purposes for fishery and survey data
# (adapted from Craig Marsh)
# Inludes explorations for: 1) Age Composition data across time, space, and sex
# 2) Length composition data across time, space, and sex
# 3) Depths that are fished across time and space (for selectivity and catchability understanding)
# 4) Regresison trees for fishery and survey lengths
# 5) Tag data explorations (Time-at-liberty distributions, time-at-liberty and distance moved, tag recovery rates, 
# tag recovoeries outside of stock boundaries, tag recoveries by region and size)
# 6) Proportion of catch relative to proportion of samples for fishery
# 7) Abundance indices and catch
# 8) Growth differences
# 9) Data comparison (quantity and time-series) by 1-area and 5-areas
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 10/2/24

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(FishFreqTree)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(geosphere)
library(SpatialSablefishAssessment)
source(here("R", 'Utility_Fxns.R'))

# Read in stat areas
nmfs_areas = read_sf(dsn = here("data", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas = st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas = nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas = st_shift_longitude(nmfs_areas) # shift longitude for plotting
nmfs_areas = nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>% 
  mutate(
    NAME = case_when(
      NAME == "Aleutian Islands" ~ "Aleutian Islands (AI)",
      NAME == "Bering Sea" ~ "Bering Sea (BS)",
      NAME == "Western Gulf of Alaska" ~ "Western Gulf of Alaska (WGOA)",
      NAME == "Central Gulf of Alaska" ~ "Central Gulf of Alaska (CGOA)",
      NAME == "West Yakutat" ~ "West Yakutat (WY)",
      NAME == "East Yakutat / Southeast Alaska" ~ "East Yakutat/Southeast (EY/SE)"
    ),
    NAME = factor(NAME, levels = c('Bering Sea (BS)', 'Aleutian Islands (AI)', 'Western Gulf of Alaska (WGOA)', 
                                   'Central Gulf of Alaska (CGOA)', 'West Yakutat (WY)', "East Yakutat/Southeast (EY/SE)")))

# Get map for plotting
west = ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west = st_shift_longitude(west) # shift ongitude for plotting
ak = subset(ne_states(c("United States of America"), returnclass = "sf"), name == 'Alaska')
ak = st_shift_longitude(ak) # shift ongitude for plotting


# Set up data labels
# region, sex, age, labels
rsa_labels = c(
  paste('BS', 'F', 2:31, sep = '-'),
  paste('BS', 'M', 2:31, sep = '-'),
  paste('AI', 'F', 2:31, sep = '-'),
  paste('AI', 'M', 2:31, sep = '-'),
  paste('WGOA', 'F', 2:31, sep = '-'),
  paste('WGOA', 'M', 2:31, sep = '-'),
  paste('CGOA', 'F', 2:31, sep = '-'),
  paste('CGOA', 'M', 2:31, sep = '-'),
  paste('WY', 'F', 2:31, sep = '-'),
  paste('WY', 'M', 2:31, sep = '-'),
  paste('EY/SE', 'F', 2:31, sep = '-'),
  paste('EY/SE', 'M', 2:31, sep = '-')
)

# region age labels
ra_labels = c(
  paste('BS', 2:31, sep = '-'),
  paste('AI', 2:31, sep = '-'),
  paste('WGOA', 2:31, sep = '-'),
  paste('CGOA', 2:31, sep = '-'),
  paste('WY', 2:31, sep = '-'),
  paste('EY/SE', 2:31, sep = '-')
)

# region, sex, length, labels
rsl_labels = c(
  paste('BS', 'F', 41:99, sep = '-'),
  paste('BS', 'M', 41:99, sep = '-'),
  paste('AI', 'F', 41:99, sep = '-'),
  paste('AI', 'M', 41:99, sep = '-'),
  paste('WGOA', 'F', 41:99, sep = '-'),
  paste('WGOA', 'M', 41:99, sep = '-'),
  paste('CGOA', 'F', 41:99, sep = '-'),
  paste('CGOA', 'M', 41:99, sep = '-'),
  paste('WY', 'F', 41:99, sep = '-'),
  paste('WY', 'M', 41:99, sep = '-'),
  paste('EY/SE', 'F', 41:99, sep = '-'),
  paste('EY/SE', 'M', 41:99, sep = '-')
)

# region length labels
rl_labels = c(
  paste('BS', 41:99, sep = '-'),
  paste('AI', 41:99, sep = '-'),
  paste('WGOA',41:99, sep = '-'),
  paste('CGOA',41:99, sep = '-'),
  paste('WY',41:99, sep = '-'),
  paste('EY/SE',41:99, sep = '-')
)

colors = unname(ggthemes::ggthemes_data[["colorblind"]][["value"]]) # get colors

# Survey ------------------------------------------------------------------
# Read in data
survey_length_df = read.csv(file.path("Data", "Survey", "length_summary_view.csv"))
survey_age_df = read.csv(file.path("Data", "Survey", "age_view.csv"))
haul_df = read.csv(file.path("Data", "Survey", "Haul.csv"))
catch_df = read.csv(file.path("Data", "Survey", "catch_summary_view.csv"))
index_by_area_df = read.csv(file.path("Data", "Survey", "RPN_index_Strat_3_7.csv"))
index_df = read.csv(file.path("Data", "Survey", "Alaska_full_area.csv"))

### Ages --------------------------------------------------------------------

# Characterize age data 
age_props_by_year_region_sex = survey_age_df %>% filter(Sex != 3, Age >= 2) %>% 
  mutate(
    Mgmt_Abrev = case_when(
      NPFMC.Sablefish.Mgmt.Area == "Aleutians" ~ "AI",
      NPFMC.Sablefish.Mgmt.Area == "Bering Sea" ~ "BS",
      NPFMC.Sablefish.Mgmt.Area == "Western Gulf of Alaska" ~ "WGOA",
      NPFMC.Sablefish.Mgmt.Area == "Central Gulf of Alaska" ~ "CGOA",
      NPFMC.Sablefish.Mgmt.Area == "West Yakutat" ~ "WY",
      NPFMC.Sablefish.Mgmt.Area == "East Yakutat/Southeast" ~ "EY/SE"
    ),
    Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', 'EY/SE')),
    Age = ifelse(Age >= 31, 31, Age),
    Sex_Abrev = ifelse(Sex == 1, "M", "F"),
    RSA_Cat = paste(Mgmt_Abrev, Sex_Abrev, Age, sep = '-'), # Region, Sex, Age
    RSA_Cat = factor(RSA_Cat, levels = rsa_labels), # Region, Sex, Age
    RA_Cat = paste(Mgmt_Abrev, Age, sep = '-'), # Region, Age
    RA_Cat = factor(RA_Cat, levels = ra_labels) # Region, Age
  ) %>% 
  drop_na()

plot_labels <- ra_labels
plot_labels[-seq(1,length(plot_labels), 10)] <- ""

# Summarize these comps
age_comp_sum = age_props_by_year_region_sex %>% 
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RA_Cat, Age) %>% 
  summarize(age_count = n()) %>% 
  group_by(Year) %>%  # normalize within year and sex
  mutate(age_comp = age_count / sum(age_count))

# By Year
# Females
ggplot(age_comp_sum %>% filter(Sex_Abrev == "F"),
       aes(x = RA_Cat, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Age Category", y = "Age Proportions", fill = "Management Region", title = 'Females')

# Males
ggplot(age_comp_sum %>% filter(Sex_Abrev == "M"),
       aes(x = RA_Cat, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Age Category", y = "Age Proportions", fill = "Management Region", title = 'Males')

# Aggregated across years
age_comp_sum = age_props_by_year_region_sex %>% 
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RA_Cat, Age) %>% 
  summarize(age_count = n()) %>% # get counts
  group_by(Mgmt_Abrev, Sex_Abrev, RA_Cat, Age) %>%
  summarize(total_age_count = sum(age_count)) %>% # sum up across these cateogries
  ungroup() %>% 
  mutate(age_comp = total_age_count/sum(total_age_count)) # normalize across ages, sexes, and regions

# Separate out by region and sex
srv_agg_ages = ggplot(age_comp_sum, aes(x = Age, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col(position = 'identity', alpha = 0.5, color = 'black') +
  facet_wrap(~Sex_Abrev) +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none', legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(Sex_Abrev~Mgmt_Abrev) +
  labs(x = "Age Category", y = "Survey Age Proportions", fill = "Management Region")


### Survey Cohort Map -------------------------------------------------------
surv_sf <- st_as_sf(survey_age_df %>% filter(!is.na(End.Longitude..DD.), !is.na(End.Latitude..DD.)), coords = c("End.Longitude..DD.", "End.Latitude..DD."), crs = 4326) # make sf
st_crs(surv_sf) <- st_crs(ak) # set same crs
srv_sf_shifted_end <- st_shift_longitude(surv_sf) # shift long
srv_long <- st_coordinates(srv_sf_shifted_end)[,1] # get adjusted long
survey_age_df$Lon2 <- srv_long

# Survey Cohort
survey_cohort <- survey_age_df %>% 
  mutate(Cohort = Year - Age) %>% 
  filter(Cohort >= 1996, Cohort < 2016) %>%
  group_by(Year, Cohort) %>% 
  drop_na() %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) %>% 
  summarize(mean_lon = mean(Lon2),
           mean_lat = mean(End.Latitude..DD.)) 

# Cohort plot
ggplot() +
geom_sf(data = ak) +
geom_text(survey_cohort, mapping = aes(x = mean_lon, y = mean_lat, label = Year, color = Year), size = 3) +
geom_path(survey_cohort, mapping = aes(x = mean_lon, y = mean_lat, color = Year, group = Cohort), size = 1, arrow = arrow()) +
scale_color_viridis_c() +
facet_wrap(~Cohort, ncol = 5) +
theme(axis.text.x = element_text(angle = 90)) +
theme_bw() +
labs(x = "Longitude", y = "Latitude", color = 'Year')


### Length ------------------------------------------------------------------

# Characterize age data
len_props_by_year_region_sex = survey_length_df %>% filter(Sex != 3, Length >= 41) %>%
  mutate(
    Mgmt_Abrev = case_when(
      NPFMC.Sablefish.Management.Area == "Aleutians" ~ "AI",
      NPFMC.Sablefish.Management.Area == "Bering Sea" ~ "BS",
      NPFMC.Sablefish.Management.Area == "Western Gulf of Alaska" ~ "WGOA",
      NPFMC.Sablefish.Management.Area == "Central Gulf of Alaska" ~ "CGOA",
      NPFMC.Sablefish.Management.Area == "West Yakutat" ~ "WY",
      NPFMC.Sablefish.Management.Area == "East Yakutat/Southeast" ~ "EY/SE"
    ),
    Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE")),
    Length = ifelse(Length >= 99, 99, Length),
    Sex_Abrev = ifelse(Sex == 1, "M", "F"),
    RSL_Cat = paste(Mgmt_Abrev, Sex_Abrev, Length, sep = '-'), # Region, Sex, Length
    RSL_Cat = factor(RSL_Cat, levels = rsl_labels), # Region, Sex, Length
    RL_Cat = paste(Mgmt_Abrev, Length, sep = '-'), # Region, Length
    RL_Cat = factor(RL_Cat, levels = rl_labels) # Region, Length
  ) %>%
  drop_na()

plot_labels <- rl_labels
plot_labels[-seq(1,length(plot_labels), 30)] <- ""

# Summarize these comps
len_comp_sum = len_props_by_year_region_sex %>%
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RL_Cat, Length) %>%
  summarize(len_count = n()) %>%
  group_by(Year) %>%  # normalize within year and sex
  mutate(len_comp = len_count / sum(len_count))

# By Year
# Females
ggplot(len_comp_sum %>% filter(Sex_Abrev == "F"),
       aes(x = RL_Cat, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Length Category", y = "Length Proportions", fill = "Management Region", title = 'Females')

# Males
ggplot(len_comp_sum %>% filter(Sex_Abrev == "M"), aes(x = RL_Cat, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Length Category", y = "Length Proportions", fill = "Management Region", title = 'Females')

# Aggregated across years
len_comp_sum = len_props_by_year_region_sex %>%
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RL_Cat, Length) %>%
  summarize(len_count = n()) %>% # get counts
  group_by(Mgmt_Abrev, Sex_Abrev, RL_Cat, Length) %>%
  summarize(total_len_count = sum(len_count)) %>% # sum up across these cateogries
  ungroup() %>%
  mutate(len_comp = total_len_count/sum(total_len_count)) # normalize across lengths, sexes, and regions

# Separate out by region and sex
ggplot(len_comp_sum, aes(x = Length, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col(position = 'identity', alpha = 1) +
  facet_wrap(~Sex_Abrev) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none', legend.background = element_blank()) +
  facet_grid(Sex_Abrev~Mgmt_Abrev) +
  labs(x = "Length Category", y = "Length Proportions", fill = "Management Region")

### Length (Regression Tree) ------------------------------------------------
# Set up directory to output stuff
dir = dir.create(here("figs", "Srv_LF_Regression_Tree"))
save_dir = here("figs", "Srv_LF_Regression_Tree") # save directory

# Munge dataset into regression tree format
survey_length_tree = survey_length_df %>% 
  dplyr::select(Sex, End.Latitude..DD., End.Longitude..DD., Length, Year) %>% 
  rename(x = End.Longitude..DD., y = End.Latitude..DD., year = Year) %>% 
  mutate(season = 'Q1', 
         length_labels = 
           case_when(
             Length < 41 ~ '0-41',
             Length < 43 ~ '41-43',
             Length < 45 ~ '43-45',
             Length < 47 ~ '45-47',
             Length < 49 ~ '47-49',
             Length < 51 ~ '49-51',
             Length < 53 ~ '51-53',
             Length < 55 ~ '53-55',
             Length < 57 ~ '55-57',
             Length < 59 ~ '57-59',
             Length < 61 ~ '59-61',
             Length < 63 ~ '61-63',
             Length < 65 ~ '63-65',
             Length < 67 ~ '65-67',
             Length < 69 ~ '67-69',
             Length < 73 ~ '71-73',
             Length < 75 ~ '73-75',
             Length < 77 ~ '75-77',
             Length < 79 ~ '77-79',
             Length < 81 ~ '79-81',
             Length < 83 ~ '81-83',
             Length < 85 ~ '83-85',
             Length < 87 ~ '85-87',
             Length < 89 ~ '87-89',
             Length < 43 ~ '41-43',
             Length < 45 ~ '43-45',
             Length < 47 ~ '45-47',
             Length < 49 ~ '47-49',
             Length < 51 ~ '49-51',
             Length < 53 ~ '51-53',
             Length < 55 ~ '53-55',
             Length < 57 ~ '55-57',
             Length < 59 ~ '57-59',
             Length < 61 ~ '59-61',
             Length < 63 ~ '61-63',
             Length < 65 ~ '63-65',
             Length < 67 ~ '65-67',
             Length < 69 ~ '67-69',
             Length < 73 ~ '71-73',
             Length < 75 ~ '73-75',
             Length < 77 ~ '75-77',
             Length < 79 ~ '77-79',
             Length < 81 ~ '79-81',
             Length < 83 ~ '81-83',
             Length < 85 ~ '83-85',
             Length < 87 ~ '85-87',
             Length < 89 ~ '87-89',
             Length < 93 ~ '91-93',
             Length < 95 ~ '93-95',
             Length < 97 ~ '95-97',
             Length < 99 ~ '97-99',
             TRUE ~ "99-99+"
           ))

# Shift longitude
surv_len_sf <- st_as_sf(survey_length_tree %>% filter(!is.na(x), !is.na(y)), coords = c("x", "y"), crs = 4326) # make sf
st_crs(surv_len_sf) <- st_crs(ak) # set same crs
srv_sf_shifted_end <- st_shift_longitude(surv_len_sf) # shift long
srv_long <- st_coordinates(srv_sf_shifted_end)[,1] # get adjusted long
survey_length_tree$x <- srv_long # replace longitude with shifted

# Do some residual cleaning up
lf_tree = survey_length_tree %>%
  mutate(x = round(x), y = round(y),
         x = as.numeric(
           paste(cut(x, breaks = seq(min(x), max(x), 2),
                     labels = seq(min(x), max(x) - 2, 2)))
         ),
         y = as.numeric(paste(
           cut(y, breaks = seq(min(y), max(y), 2),
               labels = seq(min(y), max(y) - 2, 2))
         ))) %>%   # do some quick rounding
  drop_na() %>% 
  group_by(x,y,season,year,length_labels, Sex) %>%
  count() %>%
  rename(layer = n) %>% 
  filter(Sex == 2)

lf_tree$quarter = as.numeric(substring(lf_tree$season, first = 2)) # rename quarter
lf_tree$lat = lf_tree$y # longitude
lf_tree$lon = lf_tree$x # latitude

# covert numbers to proportions
lf_tree = lf_tree %>% group_by(year, quarter, lat, lon) %>% mutate(sample_size = sum(layer, na.rm = T), length_props = layer / sample_size)
# pivot wider so its in correct format
LF_long = lf_tree %>% pivot_wider(id_cols = c(year, quarter, lat, lon, sample_size), names_from = length_labels, values_from = length_props, values_fill = 0)

# Set up dimensions
LF_res = 2 # bin width
LF_min_class = seq(from = 41, to = 97, by = LF_res) # make the last length bin aplus group
LF_max_class = seq(from = 41, to = 97, by = LF_res) # make the last length bin aplus group
LF_min_class = c(0, LF_min_class)
LF_max_class = c(LF_max_class, 200)
fcol = 6 # first column
lcol = 34 # last column
bins = LF_min_class[-1]

# Make sure things sum to 1
row_sum = apply(LF_long[,fcol:lcol],1,sum)
bad_ndx = which(abs(row_sum-1)>0.05)
length(bad_ndx)
LF_long[bad_ndx, ]
row_sum[bad_ndx]
LF_long = subset(LF_long, subset = abs(row_sum-1)<=0.05) # drop stuff that doesn't sum to 1

# Make LF map
make.lf.map(LF = LF_long, fcol, lcol, bins, save_dir = save_dir, plot_name = "make_lf_map", plot_format = "png")
make.meanl.map(LF_long, fcol, lcol, bins, save_dir, s =13)

# Run regression tree without year and quarter
LF_tree = run_regression_tree(LF = LF_long, fcol, lcol, bins, Nsplit = 3,
                              save_dir,manual = FALSE,select=NA,lat.min=1,
                              lon.min=1,year.min=1,quarter=F,year=F,
                              include_dummy=FALSE,pdf=FALSE)

LF_tree$Record
make.split.map(LF = LF_tree$LF, Nsplit = 3, save_dir = save_dir, s = 2)

# Draw lines for boundary splites
split_coords_1 = matrix(c(204, 60, 204, 40), byrow = T, ncol = 2)
split_coords_2 = matrix(c(190, 60, 190, 40), byrow = T, ncol = 2)
split_coords_3 = matrix(c(176, 60, 176, 40), byrow = T, ncol = 2)

split_lons_1 = st_cast(st_sfc(st_linestring(split_coords_1)), "LINESTRING")
split_lons_1 = st_set_crs(split_lons_1, st_crs(nmfs_areas))
split_lons_3 = st_cast(st_sfc(st_linestring(split_coords_3)), "LINESTRING")
split_lons_3 = st_set_crs(split_lons_3, st_crs(nmfs_areas))
split_lons_2 = st_cast(st_sfc(st_linestring(split_coords_2)), "LINESTRING")
split_lons_2 = st_set_crs(split_lons_2, st_crs(nmfs_areas))


## visualise these splits on a map.
world <- ne_countries(scale = "medium", returnclass = "sf")
srv_reg_tree_map = ggplot() +
  geom_sf(data = nmfs_areas, aes(fill = NAME, lty = GEN_NAME), alpha = 0.35, lwd = 0.4) +
  geom_sf(data = split_lons_1, linewidth = 1.5, lty = 2) +
  geom_sf(data = split_lons_2, linewidth = 1.5, lty = 2) +
  geom_sf(data = split_lons_3, linewidth = 1.5, lty = 2) +
  geom_sf(data = west) +
  
  # Annotate with labels
  annotate("text", label = "Alaska", x = 208, y = 65, size = 6) + # Alaska label
  annotate("text", label = "Canada", x = 224.65, y = 63, size = 6) + # Canada label
  annotate("text", label = "Russia", x = 172, y = 67, size = 6) + # Canada label
  
  coord_sf(ylim = c(44, 70.5), xlim = c(165, 235)) + # Restrict Map Area
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c("red", "blue")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5, pch = c(19, 19)))) +
  scale_linetype_discrete(guide = "none") +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.745, 0.19), legend.box = "vertical", legend.background = element_blank(),
        legend.spacing.y = unit(0.075, "cm"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "AK Sablefish Mgmt Regions", color = "Regression Tree Splits",
       title = 'Survey Lengths (Regression Tree Splits)')

  
### Survey Depth Samples ----------------------------------------------------
srv_depth_len_plot = survey_length_df %>% 
  group_by(Year, Stratum.Description, NPFMC.Sablefish.Management.Area) %>% 
  summarize(n = n()) %>% 
  mutate(NPFMC.Sablefish.Management.Area = 
           ifelse(NPFMC.Sablefish.Management.Area == "Aleutians", 'Aleutian Islands', NPFMC.Sablefish.Management.Area),
         NPFMC.Sablefish.Management.Area = factor(NPFMC.Sablefish.Management.Area,
                                                  levels = c("Bering Sea", "Aleutian Islands",
                                                             "Western Gulf of Alaska", "Central Gulf of Alaska",
                                                             "West Yakutat", "East Yakutat/Southeast"))) %>% 
  filter(!Stratum.Description %in% c("1201m +", "Unknown")) %>% 
  ggplot(aes(Stratum.Description, n, fill = NPFMC.Sablefish.Management.Area)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  labs(x = "Stratum Depth", y = "Number of Length Samples", fill = "AK Sablefish Mgmt Regions") +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9, 0.88), legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

srv_depth_len_plot

# Fishery ------------------------------------------------------------
observer_age_df = read.csv(file.path("Data", "Observer", "norpac_age_report.csv"), skip = 6) %>% 
  filter(!is.na(LonDD.End), !is.na(LatDD.End)) 
observer_length_df = read.csv(file.path("Data", "Observer", "norpac_length_report.csv"), skip = 6) %>% 
  filter(!is.na(LonDD.End), !is.na(LatDD.End)) 

### Ages --------------------------------------------------------------------

# Characterize age data 
age_props_by_year_region_sex = observer_age_df %>% filter(Sex != 3, Age >= 2) %>% 
  filter(!FMP.Subarea %in% c("PWSI", "SEI"), Sex != 'U') %>% 
  mutate(
    Mgmt_Abrev = case_when(
      FMP.Subarea == "AI" ~ "AI",
      FMP.Subarea == "BS" ~ "BS",
      FMP.Subarea == "WG" ~ "WGOA",
      FMP.Subarea == "CG" ~ "CGOA",
      FMP.Subarea == "WY" ~ "WY",
      FMP.Subarea == "SE" ~ "EY/SE"
    ),
    Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE")),
    Age = ifelse(Age >= 31, 31, Age),
    Sex_Abrev = Sex,
    RSA_Cat = paste(Mgmt_Abrev, Sex_Abrev, Age, sep = '-'), # Region, Sex, Age
    RSA_Cat = factor(RSA_Cat, levels = rsa_labels), # Region, Sex, Age
    RA_Cat = paste(Mgmt_Abrev, Age, sep = '-'), # Region, Age
    RA_Cat = factor(RA_Cat, levels = ra_labels), # Region, Age
    Gear_Abbrev = case_when(
      Gear.Description == "LONGLINER" ~ "Longliner",
      Gear.Description == "NON PELAGIC" ~ "Trawl",
      Gear.Description == "PELAGIC" ~ "Trawl",
      Gear.Description == "POT OR TRAP" ~ "Pot"
    )
  ) 

plot_labels <- ra_labels
plot_labels[-seq(1,length(plot_labels), 10)] <- ""

# Summarize these comps
age_comp_sum = age_props_by_year_region_sex %>% 
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RA_Cat, Age, Gear_Abbrev) %>% 
  summarize(age_count = n()) %>% 
  group_by(Year) %>%  # normalize within year and sex
  mutate(age_comp = age_count / sum(age_count))

# By Year
# Females
ggplot(age_comp_sum %>% filter(Sex_Abrev == "F", Gear_Abbrev == "Longliner"),
       aes(x = RA_Cat, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Age Category", y = "Age Proportions", fill = "Management Region", title = 'Females')

# Males
ggplot(age_comp_sum %>% filter(Sex_Abrev == "M", Gear_Abbrev == "Longliner"),
       aes(x = RA_Cat, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Age Category", y = "Age Proportions", fill = "Management Region", title = 'Males')

# Aggregated across years
age_comp_sum = age_props_by_year_region_sex %>% 
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RA_Cat, Age, Gear_Abbrev) %>% 
  summarize(age_count = n()) %>% # get counts
  group_by(Mgmt_Abrev, Sex_Abrev, RA_Cat, Age, Gear_Abbrev) %>%
  summarize(total_age_count = sum(age_count)) %>% # sum up across these cateogries
  ungroup() %>% 
  mutate(age_comp = total_age_count/sum(total_age_count)) # normalize across ages, sexes, and regions

# Separate out by region and sex
fish_agg_ages = ggplot(age_comp_sum %>% filter(Gear_Abbrev == "Longliner"), aes(x = Age, y = age_comp, fill = Mgmt_Abrev)) +
  geom_col(position = 'identity', alpha = 0.5, color = 'black') +
  facet_wrap(~Sex_Abrev) +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none', legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(Sex_Abrev~Mgmt_Abrev) +
  labs(x = "Age Category", y = "Fishery Age Proportions", fill = "Management Region")


### Fishery Cohort Map -------------------------------------------------------
fish_sf <- st_as_sf(observer_age_df %>% filter(!is.na(LonDD.End), !is.na(LatDD.End)), coords = c("LonDD.End", "LatDD.End"), crs = 4326) # make sf
st_crs(fish_sf) <- st_crs(ak) # set same crs
fish_sf_shifted_end <- st_shift_longitude(fish_sf) # shift long
fish_long <- st_coordinates(fish_sf_shifted_end)[,1] # get adjusted long
observer_age_df$Lon2 <- fish_long


observer_age_df %>% 
  group_by(Gear.Description) %>% count()
# Fishery Cohort
fish_cohort <- observer_age_df %>% 
  mutate(Cohort = Year - Age) %>% 
  filter(Cohort >= 1996, Cohort < 2016) %>%
  group_by(Year, Cohort) %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) %>% 
  summarize(mean_lon = mean(Lon2, na.rm = T),
            mean_lat = mean(LatDD.End, na.rm = T)) 

# Cohort plot (Longliner)
ggplot() +
  geom_text(fish_cohort, mapping = aes(x = mean_lon, y = mean_lat, label = Year, color = Year), size = 3) +
  geom_path(fish_cohort, mapping = aes(x = mean_lon, y = mean_lat, color = Year, group = Cohort), size = 1, arrow = arrow()) +
  geom_sf(data = ak) +
  scale_color_viridis_c() +
  facet_wrap(~Cohort, ncol = 5) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = 'Year')

### Length ------------------------------------------------------------------

# Characterize age data
len_props_by_year_region_sex = observer_length_df %>% filter(Sex != 3, Length..cm. >= 41) %>%
  rename(Length = Length..cm.) %>% 
  filter(!FMP.Subarea %in% c("PWSI", "SEI", "", "WOC"), Sex != 'U') %>% 
  mutate(
    Mgmt_Abrev = case_when(
      FMP.Subarea == "AI" ~ "AI",
      FMP.Subarea == "BS" ~ "BS",
      FMP.Subarea == "WG" ~ "WGOA",
      FMP.Subarea == "CG" ~ "CGOA",
      FMP.Subarea == "WY" ~ "EGOA",
      FMP.Subarea == "SE" ~ "EGOA"
    ),
    Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
    Length = ifelse(Length >= 99, 99, Length),
    Sex_Abrev = Sex,
    RSL_Cat = paste(Mgmt_Abrev, Sex_Abrev, Length, sep = '-'), # Region, Sex, Length
    RSL_Cat = factor(RSL_Cat, levels = rsl_labels), # Region, Sex, Length
    RL_Cat = paste(Mgmt_Abrev, Length, sep = '-'), # Region, Length
    RL_Cat = factor(RL_Cat, levels = rl_labels), # Region, Length
    Gear_Abbrev = ifelse(Gear.Description == "LONGLINER", 'Longliner', 
                         ifelse(Gear.Description == 'POT OR TRAP', "Pot", "Trawl"))
  ) 

plot_labels <- rl_labels
plot_labels[-seq(1,length(plot_labels), 30)] <- ""

# Summarize these comps
len_comp_sum = len_props_by_year_region_sex %>%
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RL_Cat, Length, Gear_Abbrev) %>%
  summarize(len_count = n()) %>%
  group_by(Year) %>%  # normalize within year and sex
  mutate(len_comp = len_count / sum(len_count))

# By Year
# Females
ggplot(len_comp_sum %>% filter(Sex_Abrev == "F", Gear_Abbrev == "Longliner"),
       aes(x = RL_Cat, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Length Category", y = "Length Proportions", fill = "Management Region", title = 'Females')

# Males
ggplot(len_comp_sum %>% filter(Sex_Abrev == "M", Gear_Abbrev == "Longliner"), 
       aes(x = RL_Cat, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col() +
  facet_wrap(~Sex_Abrev) +
  scale_x_discrete(labels = plot_labels) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top", legend.background = element_blank()) +
  facet_wrap(~Year) +
  labs(x = "Region, Length Category", y = "Length Proportions", fill = "Management Region", title = 'Females')

# Aggregated across years
len_comp_sum = len_props_by_year_region_sex %>%
  group_by(Year, Mgmt_Abrev, Sex_Abrev, RL_Cat, Length, Gear_Abbrev) %>%
  summarize(len_count = n()) %>% # get counts
  group_by(Mgmt_Abrev, Sex_Abrev, RL_Cat, Length, Gear_Abbrev) %>%
  summarize(total_len_count = sum(len_count)) %>% # sum up across these cateogries
  ungroup() %>%
  mutate(len_comp = total_len_count/sum(total_len_count)) # normalize across lengths, sexes, and regions

# Separate out by region and sex
ggplot(len_comp_sum %>% filter(Gear_Abbrev == "Longliner"), aes(x = Length, y = len_comp, fill = Mgmt_Abrev)) +
  geom_col(position = 'identity', alpha = 1) +
  facet_wrap(~Sex_Abrev) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none', legend.background = element_blank()) +
  facet_grid(Sex_Abrev~Mgmt_Abrev) +
  labs(x = "Length Category", y = "Length Proportions", fill = "Management Region")


### Length (Regression Tree) ------------------------------------------------

# Set up directory to output stuff
dir = dir.create(here("figs", "Fish_LF_Regression_Tree"))
save_dir = here("figs", "Fish_LF_Regression_Tree") # save directory

# Munge dataset into regression tree format
fish_length_tree = observer_length_df %>% 
  mutate(Gear.Description = ifelse(Gear.Description %in% c("LONGLINER", "POT OR TRAP"), "Fixed", "Trawl")) %>% 
  filter(Gear.Description == "Fixed") %>% 
  select(Sex, LatDD.End, LonDD.End, Length..cm., Year) %>% 
  rename(x = LonDD.End, y = LatDD.End, year = Year, Length = Length..cm.) %>% 
  mutate(season = 'Q1', 
         length_labels = 
           case_when(
             Length < 41 ~ '0-41',
             Length < 43 ~ '41-43',
             Length < 45 ~ '43-45',
             Length < 47 ~ '45-47',
             Length < 49 ~ '47-49',
             Length < 51 ~ '49-51',
             Length < 53 ~ '51-53',
             Length < 55 ~ '53-55',
             Length < 57 ~ '55-57',
             Length < 59 ~ '57-59',
             Length < 61 ~ '59-61',
             Length < 63 ~ '61-63',
             Length < 65 ~ '63-65',
             Length < 67 ~ '65-67',
             Length < 69 ~ '67-69',
             Length < 73 ~ '71-73',
             Length < 75 ~ '73-75',
             Length < 77 ~ '75-77',
             Length < 79 ~ '77-79',
             Length < 81 ~ '79-81',
             Length < 83 ~ '81-83',
             Length < 85 ~ '83-85',
             Length < 87 ~ '85-87',
             Length < 89 ~ '87-89',
             Length < 43 ~ '41-43',
             Length < 45 ~ '43-45',
             Length < 47 ~ '45-47',
             Length < 49 ~ '47-49',
             Length < 51 ~ '49-51',
             Length < 53 ~ '51-53',
             Length < 55 ~ '53-55',
             Length < 57 ~ '55-57',
             Length < 59 ~ '57-59',
             Length < 61 ~ '59-61',
             Length < 63 ~ '61-63',
             Length < 65 ~ '63-65',
             Length < 67 ~ '65-67',
             Length < 69 ~ '67-69',
             Length < 73 ~ '71-73',
             Length < 75 ~ '73-75',
             Length < 77 ~ '75-77',
             Length < 79 ~ '77-79',
             Length < 81 ~ '79-81',
             Length < 83 ~ '81-83',
             Length < 85 ~ '83-85',
             Length < 87 ~ '85-87',
             Length < 89 ~ '87-89',
             Length < 93 ~ '91-93',
             Length < 95 ~ '93-95',
             Length < 97 ~ '95-97',
             Length < 99 ~ '97-99',
             TRUE ~ "99-99+"
           ))

# Shift longitude
fish_len_sf <- st_as_sf(fish_length_tree %>% filter(!is.na(x), !is.na(y)), coords = c("x", "y"), crs = 4326) # make sf
st_crs(fish_len_sf) <- st_crs(ak) # set same crs
fish_sf_shifted_end <- st_shift_longitude(fish_len_sf) # shift long
fish_long <- st_coordinates(fish_sf_shifted_end)[,1] # get adjusted long
fish_length_tree$x <- fish_long # replace longitude with shifted

# Do some residual cleaning up
lf_tree = fish_length_tree %>%
  mutate(x = round(x), y = round(y),
         x = as.numeric(
           paste(cut(x, breaks = seq(min(x), max(x), 2),
                     labels = seq(min(x), max(x) - 2, 2)))
         ),
         y = as.numeric(paste(
           cut(y, breaks = seq(min(y), max(y), 2),
               labels = seq(min(y), max(y) - 2, 2))
         ))) %>%   # do some quick rounding
  group_by(x,y,season,year,length_labels,Sex) %>%
  count() %>%
  rename(layer = n) %>% 
  filter(Sex == 'M') %>% 
  drop_na() %>% 
  distinct()

lf_tree$quarter = as.numeric(substring(lf_tree$season, first = 2)) # rename quarter
lf_tree$lat = lf_tree$y # longitude
lf_tree$lon = lf_tree$x # latitude

# covert numbers to proportions
lf_tree = lf_tree %>% 
  group_by(year, quarter, lat, lon) %>%
  mutate(sample_size = sum(layer, na.rm = T), length_props = layer / sample_size)

# pivot wider so its in correct format
LF_long = lf_tree %>% pivot_wider(id_cols = c(year, quarter, lat, lon, sample_size), 
                                  names_from = length_labels, values_from = length_props, values_fill = 0)

# Set up dimensions
LF_res = 2 # bin width
LF_min_class = seq(from = 41, to = 97, by = LF_res) # make the last length bin aplus group
LF_max_class = seq(from = 41, to = 97, by = LF_res) # make the last length bin aplus group
LF_min_class = c(0, LF_min_class)
LF_max_class = c(LF_max_class, 200)
fcol = 6 # first column
lcol = 34 # last column
bins = LF_min_class[-1]

# Make sure things sum to 1
row_sum = apply(LF_long[,fcol:lcol],1,sum)
bad_ndx = which(abs(row_sum-1)>0.05)
length(bad_ndx)
LF_long[bad_ndx, ]
row_sum[bad_ndx]
LF_long = subset(LF_long, subset = abs(row_sum-1)<=0.05) # drop stuff that doesn't sum to 1

# Make LF map
make.lf.map(LF = LF_long, fcol, lcol, bins, save_dir = save_dir, plot_name = "make_lf_map", plot_format = "png")
make.meanl.map(LF_long, fcol, lcol, bins, save_dir, s =13)

# Run regression tree without year and quarter
LF_tree = run_regression_tree(LF = LF_long, fcol, lcol, bins, Nsplit = 3,
                              save_dir,manual = FALSE,select=NA,lat.min=1,
                              lon.min=1,year.min=1,quarter=F,year=F,
                              include_dummy=FALSE,pdf=FALSE)

LF_tree$Record
make.split.map(LF = LF_tree$LF, Nsplit = 3, save_dir = save_dir, s = 2)

# Draw lines for boundary splites
split_coords_1 = matrix(c(180, 60, 180, 40), byrow = T, ncol = 2)
split_coords_2 = matrix(c(212, 60, 212, 40), byrow = T, ncol = 2)
split_coords_3 = matrix(c(178, 60, 178, 40), byrow = T, ncol = 2)

split_lons_1 = st_cast(st_sfc(st_linestring(split_coords_1)), "LINESTRING")
split_lons_1 = st_set_crs(split_lons_1, st_crs(nmfs_areas))
split_lons_3 = st_cast(st_sfc(st_linestring(split_coords_3)), "LINESTRING")
split_lons_3 = st_set_crs(split_lons_3, st_crs(nmfs_areas))
split_lons_2 = st_cast(st_sfc(st_linestring(split_coords_2)), "LINESTRING")
split_lons_2 = st_set_crs(split_lons_2, st_crs(nmfs_areas))


## visualise these splits on a map.
world <- ne_countries(scale = "medium", returnclass = "sf")
fsh_reg_tree_map = ggplot() +
  geom_sf(data = nmfs_areas, aes(fill = NAME, lty = GEN_NAME), alpha = 0.35, lwd = 0.4) +
  geom_sf(data = split_lons_1, linewidth = 1.5, lty = 2) +
  geom_sf(data = split_lons_2, linewidth = 1.5, lty = 2) +
  geom_sf(data = split_lons_3, linewidth = 1.5, lty = 2) +
  geom_sf(data = west) +
  scale_linetype_discrete(guide = "none") +
  scale_fill_manual(values = colors) +
  # Annotate with labels
  annotate("text", label = "Alaska", x = 208, y = 65, size = 6) + # Alaska label
  annotate("text", label = "Canada", x = 224.65, y = 63, size = 6) + # Canada label
  annotate("text", label = "Russia", x = 172, y = 67, size = 6) + # Canada label
  coord_sf(ylim = c(44, 70.5), xlim = c(165, 235)) + # Restrict Map Area
  theme_bw(base_size = 13) +
  theme(legend.position = 'none', 
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Longitude", y = "Latitude", fill = "AK Sablefish Mgmt Regions", color = "Regression Tree Splits",
       title = 'Fishery Lengths (Regression Tree Splits)')

### Fishery Depth and Sampling Locations ---------------------------------------------------
obs_plot_df = observer_length_df %>% 
  mutate(Gear.Description = ifelse(!Gear.Description %in% c("LONGLINER", "POT OR TRAP"), "TRAWL", Gear.Description)) %>% 
  filter(!FMP.Subarea %in% c("PWSI", "SEI"), Sex != 'U') %>% 
  mutate(
    Mgmt_Abrev = case_when(
      FMP.Subarea == "AI" ~ "AI",
      FMP.Subarea == "BS" ~ "BS",
      FMP.Subarea == "WG" ~ "WGOA",
      FMP.Subarea == "CG" ~ "CGOA",
      FMP.Subarea == "WY" ~ "WY",
      FMP.Subarea == "SE" ~ "EY/SE"
    ),
    Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE"))) %>% 
  filter(!is.na(Mgmt_Abrev))

# Depth
bot_depth_len_obs = ggplot(obs_plot_df, aes(Year, Bottom.Depth..Fathoms. * 1.829, group = Year)) +
geom_boxplot(outlier.colour = NA, alpha = 0.5, outliers = F) +
geom_hline(obs_plot_df %>% group_by(Gear.Description, Mgmt_Abrev) %>% 
           summarize(Bottom.Depth..Fathoms. = mean((Bottom.Depth..Fathoms. * 1.829), na.rm = T)),
           mapping = aes(yintercept = Bottom.Depth..Fathoms.), lty = 2, color = 'blue', lwd = 1.3) +
labs(x = "Year", y = "Bottom Depth (m) of Fishing Events") +
facet_grid(Gear.Description~Mgmt_Abrev) +
ggthemes::scale_fill_colorblind() +
theme_bw(base_size = 15) +
theme(plot.background = element_rect(fill = "transparent", colour = NA))

bot_depth_len_obs

# End Longitude of fishing locations
fish_sf <- st_as_sf(obs_plot_df %>% filter(!is.na(LonDD.End), !is.na(LatDD.End)), coords = c("LonDD.End", "LatDD.End"), crs = 4326) # make sf
st_crs(fish_sf) <- st_crs(ak) # set same crs
fish_sf_shifted_end <- st_shift_longitude(fish_sf) # shift long
fish_long <- st_coordinates(fish_sf_shifted_end)[,1] # get adjusted long
obs_plot_df$x <- fish_long # replace longitude with shifted

# Longtiude 
lon_len_obs = ggplot(obs_plot_df, aes(Year, x, group = Year)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, outliers = F) +
  geom_hline(obs_plot_df %>% group_by(Gear.Description, Mgmt_Abrev) %>% 
               summarize(x = mean(x, na.rm = T)),
             mapping = aes(yintercept = x), lty = 2, color = 'blue', lwd = 1.3) +
  labs(x = "Year", y = "End Longitude (Shifted) of Fishing Events") +
  facet_wrap(Gear.Description~Mgmt_Abrev, nrow = 3, scales = 'free_y') +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 15) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

lon_len_obs

# Latitude
lat_len_obs = ggplot(obs_plot_df, aes(Year, LatDD.End, group = Year)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, outliers = F) +
  geom_hline(obs_plot_df %>% group_by(Gear.Description, Mgmt_Abrev) %>% 
               summarize(LatDD.End = mean(LatDD.End, na.rm = T)),
             mapping = aes(yintercept = LatDD.End), lty = 2, color = 'blue', lwd = 1.3) +
  labs(x = "Year", y = "End Latitude of Fishing Events") +
  facet_wrap(Gear.Description~Mgmt_Abrev, nrow = 3, scales = 'free_y') +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 15) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

lat_len_obs

# Tagging Data ------------------------------------------------------------
haul_df = read.csv(file = file.path(here("data", "Tagging"), "HAUL_VIEW.csv"))
release_df = read.csv(file = file.path(here("data", "Tagging"), "RELEASE_MVIEW_Sablefish.csv"))
recovery_df = read.csv(file = file.path(here("data", "Tagging"), "RECOVERY_MVIEW_Sablefish.csv"))

# Read in stat areas
nmfs_areas = read_sf(dsn = here("data", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas = st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas = nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas = st_shift_longitude(nmfs_areas) # shift longitude for plotting
nmfs_areas = nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>% 
  mutate(
    NAME = case_when(
      NAME == "Aleutian Islands" ~ "Aleutian Islands (AI)",
      NAME == "Bering Sea" ~ "Bering Sea (BS)",
      NAME == "Western Gulf of Alaska" ~ "Western Gulf of Alaska (WGOA)",
      NAME == "Central Gulf of Alaska" ~ "Central Gulf of Alaska (CGOA)",
      NAME == "West Yakutat" ~ "West Yakutat (EGOA)",
      NAME == "East Yakutat / Southeast Alaska" ~ "East Yakutat/Southeast (EGOA)"
    ),
    NAME = factor(NAME, levels = c('Bering Sea (BS)', 'Aleutian Islands (AI)', 'Western Gulf of Alaska (WGOA)', 
                                   'Central Gulf of Alaska (CGOA)', 'West Yakutat (EGOA)', "East Yakutat/Southeast (EGOA)")))
# Get map for plotting
west = ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west = st_shift_longitude(west) # shift ongitude for plotting

### Recoveries ----------------------------------------------------------
# shift longitude
recovery_sf_start <- st_as_sf(recovery_df %>% filter(!is.na(HLNG), !is.na(HLAT)), coords = c("HLNG", "HLAT"), crs = 4326) # make sf 
recovery_sf_end <- st_as_sf(recovery_df %>% filter(!is.na(RLNG), !is.na(RLAT)), coords = c("RLNG", "RLAT"), crs = 4326) # make sf
st_crs(recovery_sf_start) <- st_crs(nmfs_areas) # set same crs
st_crs(recovery_sf_end) <- st_crs(nmfs_areas) # set same crs
recovery_sf_shifted_start <- st_shift_longitude(recovery_sf_start) # shift long
recovery_sf_shifted_end <- st_shift_longitude(recovery_sf_end) # shift long
shift_lon_start <- st_coordinates(recovery_sf_shifted_start)[,1] # get adjusted long
shift_lon_end <- st_coordinates(recovery_sf_shifted_end)[,1] # get adjust long
recovery_df$HLNG2 <- shift_lon_start
recovery_df$RLNG2 <- shift_lon_end

# Remove points that are not in management boundaries
# inside_boundaries = st_intersects(recovery_sf_end, nmfs_areas, sparse = FALSE) # find points that are not inside any boundaries
# recovery_df = recovery_df[rowSums(inside_boundaries) == 1, ] # filter to those that are within a boundary

# Define east and west general direction
recovery_df = recovery_df %>% 
  filter(TIME_OUT > 0, RLAT != 0) %>% 
  mutate(Direction = ifelse(RLNG2 > HLNG2, "East", "West"))

# Get centroids
mean_sd_centroid = recovery_df %>% 
  summarize(mean_HLNG = mean(HLNG2),
            SD_HLNG = sd(HLNG2),
            mean_HLAT = mean(HLAT),
            SD_HLAT = sd(HLAT),
            mean_RLNG = mean(RLNG2),
            SD_RLNG = sd(RLNG2),
            mean_RLAT = mean(RLAT),
            SD_RLAT = sd(RLAT))

recovery_map_plot = ggplot() +
  geom_sf(data = nmfs_areas, aes(fill = NAME, lty = GEN_NAME), alpha = 0.35, lwd = 0.2) + # Sablefish Stat Area
  geom_segment(recovery_df, mapping = aes(x = HLNG2, HLAT, xend = RLNG2, yend = RLAT), alpha = 0.015, lineend = "round") + # Direction
  geom_point(recovery_df, mapping = aes(x = HLNG2, y = HLAT), pch = 20, position = position_jitter(width = 0.05), alpha = 0.1, size = 2) + # Release
  geom_point(recovery_df, mapping = aes(x = RLNG2, y = RLAT, color = Direction), size = 2, alpha = 0.1, pch = 18, position = position_jitter(width = 0.05)) + # Recovery 
  geom_point(mean_sd_centroid, mapping = aes(x = mean_HLNG, y = mean_HLAT), pch = 19, size = 5, color = 'white') + # Release centroid
  geom_errorbar(mean_sd_centroid, mapping = aes(x = mean_HLNG, y = mean_HLAT,  # Release SD
                                                ymin = mean_HLAT - SD_HLAT, ymax = mean_HLAT + SD_HLAT), lty = 2, color = 'white', lwd = 0.55) +
  geom_errorbarh(mean_sd_centroid, mapping = aes(y = mean_HLAT, xmin = mean_HLNG - SD_HLNG, # Release SD
                                                 xmax = mean_HLNG + SD_HLNG), lty = 2, color = 'white', lwd = 0.55) +
  geom_sf(data = west, lwd = 0.2, color = 'black') + # World Map
  
  # Annotate with labels
  annotate("text", label = "Alaska", x = 208, y = 65, size = 8) + # Alaska label
  annotate("text", label = "Canada", x = 224.65, y = 63, size = 8) + # Canada label
  annotate("text", label = "Russia", x = 172, y = 67, size = 8) + # Canada label
  
  coord_sf(ylim = c(44, 70.5), xlim = c(165, 235)) + # Restrict Map Area
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c("red", "blue")) +
  scale_linetype_discrete(guide = "none") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5, pch = c(19, 19)))) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.12, 0.155), legend.box = "vertical", legend.background = element_blank(),
        legend.spacing.y = unit(0.075, "cm"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = "Longitude", y = "Latitude", color = "Tag Recovery Direction", fill = "AK Sablefish Mgmt Regions")

### Time at Liberty + Tag Recoveries + Recovery Rates ----------------------------------------
haul_df = read.csv(file = file.path(here("data", "Tagging"), "HAUL_VIEW.csv"))
release_df = read.csv(file = file.path(here("data", "Tagging"), "RELEASE_MVIEW_Sablefish.csv"))
recovery_df = read.csv(file = file.path("Data", "Tagging", "RECOVERY_MVIEW_Sablefish_Groomed for Craig.csv"))

# shift longitude
recovery_sf_start <- st_as_sf(recovery_df %>% filter(!is.na(HLNG), !is.na(HLAT)), coords = c("HLNG", "HLAT"), crs = 4326) # make sf 
recovery_sf_end <- st_as_sf(recovery_df %>% filter(!is.na(RLNG), !is.na(RLAT)), coords = c("RLNG", "RLAT"), crs = 4326) # make sf
st_crs(recovery_sf_start) <- st_crs(nmfs_areas) # set same crs
st_crs(recovery_sf_end) <- st_crs(nmfs_areas) # set same crs
recovery_sf_shifted_start <- st_shift_longitude(recovery_sf_start) # shift long
recovery_sf_shifted_end <- st_shift_longitude(recovery_sf_end) # shift long
shift_lon_start <- st_coordinates(recovery_sf_shifted_start)[,1] # get adjusted long
shift_lon_end <- st_coordinates(recovery_sf_shifted_end)[,1] # get adjust long
recovery_df$HLNG2 <- shift_lon_start
recovery_df$RLNG2 <- shift_lon_end

# Remove points that are not in management boundaries
inside_boundaries = st_intersects(recovery_sf_end, nmfs_areas, sparse = FALSE) # find points that are not inside any boundaries
recovery_df = recovery_df[rowSums(inside_boundaries) == 1, ] # filter to those that are within a boundary

# Gear fields and other munging
recovery_plot_df = recovery_df %>% 
  filter(TIME_OUT > 0, !is.na(HLNG), 
         !is.na(HLAT), !is.na(RLNG), !is.na(RLAT)) %>% # remove points wiht no start and end long lat
  mutate(GEAR = case_when(
    GEAR == 901 ~ "Pot",
    GEAR == 902 ~ "Hook-and-Line",
    TRUE ~ "Other"
  ), # Munge gear fields
  HAUL_DATE = dmy(HAUL_DATE),
  REC_DATE = dmy(REC_DATE),
  Time_Liberty = (REC_DATE - HAUL_DATE) / 365,
  Time_Liberty_Block = case_when(
    round(Time_Liberty) < 6 ~ "0 - 5 Years",
    round(Time_Liberty) < 11 ~ "6 - 10 Years",
    round(Time_Liberty) < 16 ~ "11 - 15 Years",
    round(Time_Liberty) < 21 ~ "16 - 20 Years",
    round(Time_Liberty) < 26 ~ "21 - 25 Years",
    round(Time_Liberty) < 31 ~ "26 - 30 Years",
    round(Time_Liberty) < 36 ~ "31 - 35 Years",
    TRUE ~ "36 - 40+ Years"
  ), 
  Time_Liberty_Block = factor(Time_Liberty_Block, levels = c("0 - 5 Years", 
                                                             "6 - 10 Years",
                                                             "11 - 15 Years",
                                                             "16 - 20 Years",
                                                             "21 - 25 Years",
                                                             "26 - 30 Years",
                                                             "31 - 35 Years", "36 - 40+ Years")), 
  Direction = ifelse(RLNG2 > HLNG2, "East", "West"))

# Summarize time series
tag_recovery_ts = recovery_plot_df %>% 
  group_by(REC_YEAR, GEAR) %>% 
  summarize(n = n()) %>% 
  mutate(Fishery_Reg = case_when(
    REC_YEAR < 1995 ~ NA,
    REC_YEAR < 2016 ~ 1995,
    TRUE ~ 2017
  )) %>% 
  ggplot(aes(x = REC_YEAR, y = n, color = GEAR, lty = GEAR)) +
  geom_line(lwd = 1) +
  geom_vline(mapping = aes(xintercept = Fishery_Reg), lty = 2, lwd = 0.85) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.12, 0.88), 
        legend.background = element_blank(), legend.box = element_blank(),
        legend.spacing.x = unit(1, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Recovery Year", y = "Number of Recoveries", color = "Gear Type", lty = "Gear Type")

tag_recovery_ts

# Time at liberty
liberty_hist = ggplot(recovery_plot_df, aes(x = Time_Liberty)) +
  geom_histogram(bins = 40, color = 'black', alpha = 0.5) +
  theme_bw(base_size = 13) + 
  theme(legend.position = c(0.7, 0.86), 
        legend.background = element_blank(), legend.box = element_blank(),
        legend.spacing.x = unit(1, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
  labs(x = "Time at Liberty (Years)", y = "Number of Recoveries")

liberty_hist

# Time at liberty broken up into blocks (map plot)
liberty_map_plot = ggplot(recovery_plot_df) +
  geom_sf(data = nmfs_areas, aes(fill = NAME, lty = GEN_NAME), alpha = 0.5, lwd = 0.2) + # Sablefish Stat Area
  geom_segment(recovery_plot_df, mapping = aes(x = HLNG2, HLAT, xend = RLNG2, yend = RLAT), alpha = 0.1, lineend = "round") + # Direction
  geom_point(recovery_plot_df, mapping = aes(x = HLNG2, y = HLAT), pch = 20, position = position_jitter(width = 0.05), alpha = 0.3, size = 2) + # Release
  geom_point(recovery_plot_df, mapping = aes(x = RLNG2, y = RLAT, color = Direction), size = 2, alpha = 0.3, pch = 18, position = position_jitter(width = 0.05)) + # Recovery
  geom_sf(data = west, lwd = 0.2, color = 'black') + # World Map
  coord_sf(ylim = c(48, 70.5), xlim = c(170, 227)) + # Restrict Map Area
  facet_wrap(~Time_Liberty_Block, ncol = 2) +
  # Annotate with labels
  annotate("text", label = "Alaska", x = 208, y = 65, size = 2) + # Alaska label
  annotate("text", label = "Canada", x = 224.65, y = 63, size = 2) + # Canada label
  annotate("text", label = "Russia", x = 172.5, y = 67, size = 2) + # Canada label
  scale_fill_manual(values = c("black", "#56B4E9", "#009E73","#F0E442", "#D55E00", "#E69F00")) +
  scale_linetype_discrete(guide = "none") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5, pch = c(19, 19)))) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none') +
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "Longitude", y = "Latitude", fill = "AK Sablefish Management Regions", color = "Tag Recovery Direction")

liberty_map_plot

### Distribution of Year at Liberty -----------------------------------------
recovery_df = read.csv(file = file.path("Data", "Tagging", "RECOVERY_MVIEW_Sablefish_Groomed for Craig.csv"))
release_df = read.csv(file = file.path("Data", "Tagging", "RELEASE_MVIEW_Sablefish.csv"))
haul_df = data.table::fread(file.path("Data", "Tagging", "HAUL_VIEW.csv"))

##### Residual Munging --------------------------------------------------------
release_df$abundance = 1
recovery_df$abundance = 1
## start grooming records now
start_records = nrow(release_df)
## these are records to save the effects of grooming rules
start_ndx = rep(T, nrow(release_df))
release_grooming_record = record_grooming_rule(df = release_df, index = start_ndx, catch.col = "abundance", record = NULL, rule = "Init", attribute  = "events") 
## start grooming records now
start_records = nrow(recovery_df)
## these are records to save the effects of grooming rules
start_ndx = rep(T, nrow(recovery_df))
recovery_grooming_record = record_grooming_rule(df = recovery_df, index = start_ndx, catch.col = "abundance", record = NULL, rule = "Init", attribute  = "events") 

## observed data projection decimal degrees, -180 -180 
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
km_proj = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

## link the haul and release id  Vessel, Cruise, and Haul fields
haul_df$release_id = paste0(haul_df$VESSEL, "-", haul_df$CRUISE, "-", haul_df$HAUL)
release_df$release_id = paste0(release_df$VESSEL, "-", release_df$CRUISE, "-", release_df$HAUL)
## link the recover and release id  TagType , TagNum,
release_df$recovery_id = paste0(release_df$TAG_TYPE, "-", release_df$TAG_NUM)
recovery_df$recovery_id = paste0(recovery_df$TAG_TYPE, "-", recovery_df$TAG_NUM)

## do some simple grooming
# 1) drop data with lat or long as NA
pre_count = nrow(haul_df)
haul_df = haul_df %>% filter(!is.na(HLNG), !is.na(HLAT))
cat("dropped ", pre_count - nrow(haul_df), " records \n")

dup_rel_ndx = !is.na(release_df$HLAT) & !is.na(release_df$HLNG) 
release_grooming_record = record_grooming_rule(df = release_df, index = dup_rel_ndx, catch.col = "abundance", record = release_grooming_record, rule = "NA lat and long", attribute  = "events") 
release_df = apply_grooming_rule(release_df, index = dup_rel_ndx)

dup_ndx = !is.na(recovery_df$RLNG) & !is.na(recovery_df$RLAT) 
recovery_grooming_record = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance", record = recovery_grooming_record, rule = "NA lat and long", attribute  = "events") 
recovery_df = apply_grooming_rule(recovery_df, dup_ndx)


## are all recovery fish in release fish
table(recovery_df$recovery_id %in% release_df$recovery_id)
## reformat date
haul_df$release_date = as.Date(haul_df$HAUL_DATE, format = "%d-%b-%y")
haul_df$release_year = format(haul_df$release_date, format = "%Y")
rel_years = sort(unique(haul_df$release_year))

recovery_df$recovery_date = as.Date(recovery_df$REC_DATE, format = "%d-%b-%y")
recovery_df$recovery_year = format(recovery_df$recovery_date, format = "%Y")

recovery_df$survey_recovery = !is.na(recovery_df$RCRUISE) & recovery_df$RCRUISE != "200"
print(recovery_df %>% group_by(recovery_year) %>% summarise(sum(survey_recovery)), n = 51)
recovery_df$recovery_lon = recovery_df$RLNG
recovery_df$recovery_lat = recovery_df$RLAT

## combine haul and release data frame
release_df = release_df %>% left_join(haul_df, by = "release_id", keep = F)
release_df$release_lat = release_df$HLAT.x
release_df$release_lon = release_df$HLNG.x


## start grooming records by year at this point
start_records = nrow(release_df)
## these are records to save the effects of grooming rules
start_ndx = rep(T, nrow(release_df))
release_grooming_record_rel_yr = record_grooming_rule(df = release_df, index = !start_ndx, catch.col = "abundance", record = NULL, rule = "Init", attribute  = "events", year.col = "release_year") 
## start grooming records now
start_records = nrow(recovery_df)
## these are records to save the effects of grooming rules
start_ndx = rep(T, nrow(recovery_df))
recovery_grooming_record_rec_yr = record_grooming_rule(df = recovery_df, index = !start_ndx, catch.col = "abundance", record = NULL, rule = "Init", attribute  = "events", year.col = "recovery_year") 

## combine release data frame with recovery data based on recovery id
table(duplicated(recovery_df$recovery_id)) ## duplicated recovery id's BEFORE merging release infor
recovery_df = recovery_df %>% left_join(release_df, by = "recovery_id", keep = F)
## this adds 1030 records. THis is because there are duplicated recovery_id why perhaps something to do with duplicates
table(duplicated(recovery_df$recovery_id)) ## duplicated recovery id's BEFORE merging release infor

dup_rel_ndx = !duplicated(release_df$recovery_id)
table((release_df %>% filter(!dup_rel_ndx))$HAUL_YEAR.x)

release_grooming_record = record_grooming_rule(df = release_df, index = dup_rel_ndx, catch.col = "abundance", record = release_grooming_record, rule = "Duplicated recovery id", attribute  = "events") 
release_grooming_record_rel_yr = record_grooming_rule(df = release_df, index = dup_rel_ndx, catch.col = "abundance", record = release_grooming_record_rel_yr, rule = "Duplicated recovery id", attribute  = "events", year.col = "release_year") 
release_df = apply_grooming_rule(release_df, index = dup_rel_ndx)

dup_ndx = !duplicated(recovery_df$recovery_id)
table((recovery_df %>% filter(!dup_ndx))$release_year)
recovery_grooming_record = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance.x", record = recovery_grooming_record, rule = "Duplicated recovery id", attribute  = "events") 
recovery_grooming_record_rec_yr = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance.x", record = recovery_grooming_record_rec_yr, rule = "Duplicated recovery id", attribute  = "events", year.col = "recovery_year")  
recovery_df = apply_grooming_rule(recovery_df, dup_ndx)

## release lat lon NAs in recovery data (probably should cut these out before we merge)
dup_ndx = !is.na(recovery_df$release_lat) & !is.na(recovery_df$release_lon) 
recovery_grooming_record = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance.x", record = recovery_grooming_record, rule = "NA release lat and long", attribute  = "events") 
recovery_df = apply_grooming_rule(recovery_df, dup_ndx)


## calculate time at liberty
recovery_df$time_at_liberty = recovery_df$recovery_date - recovery_df$release_date 
## calculate change in growth
recovery_df$growth_at_liberty = recovery_df$RSIZE - recovery_df$HSIZE.y
head(cbind(recovery_df$GROWTH, recovery_df$growth_at_liberty, recovery_df$HSIZE.y))
## NOTE: Becareful not to use GROWTH!! it has been calculated assuming NA = 0
## How many recapture length records do we have (as a percentage)
sum(!is.na(recovery_df$growth_at_liberty)) / nrow(recovery_df ) * 100
## quick plot of growth
hist(as.numeric(recovery_df$time_at_liberty))
summary(as.numeric(recovery_df$time_at_liberty))

## drop negative time at liberty
dup_ndx = recovery_df$time_at_liberty > 0
recovery_grooming_record = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance.x", record = recovery_grooming_record, rule = "negative time-at-liberty recovery id", attribute  = "events") 
recovery_grooming_record_rec_yr = record_grooming_rule(df = recovery_df, index = dup_ndx, catch.col = "abundance.x", record = recovery_grooming_record_rec_yr, rule = "negative time-at-liberty recovery id", attribute  = "events", year.col = "recovery_year")  
recovery_df = apply_grooming_rule(recovery_df, dup_ndx)

# now as a function of years at liberty
recovery_df$years_at_liberty = as.numeric(recovery_df$recovery_year) - as.numeric(recovery_df$release_year)
recovery_df$years_at_liberty_num = recovery_df$years_at_liberty
recovery_df = recovery_df %>% mutate(years_at_liberty = ifelse(years_at_liberty >= 15, "15+", years_at_liberty), years_at_liberty_num = ifelse(years_at_liberty_num >= 15, 15, years_at_liberty_num))
recovery_df$years_at_liberty = paste0("Year at Liberty = ", recovery_df$years_at_liberty)
recovery_df$years_at_liberty = factor(recovery_df$years_at_liberty, levels = paste0("Year at Liberty = ", c(0:15, "15+")))

##### Plot --------------------------------------------------------------------
dist = ggplot() + 
  geom_col(data = recovery_df %>% st_drop_geometry() %>%
             group_by(release_year, years_at_liberty) %>%
             summarise(recoveries = n()) %>%
             group_by(years_at_liberty) %>%
             mutate(recoveries = recoveries / max(recoveries)),
           aes(x = release_year, y = recoveries), alpha = 0.55, color = 'black') + # relative number of recoveries for a given year at liberty
  labs(y = "Distribution of Relative Number of Recoveries", x = "Release Year", col = "", fill = "") +
  facet_wrap(~years_at_liberty, ncol = 2) +
  scale_x_discrete(breaks = every_nth(n = 6)) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none', 
        legend.background = element_blank(), legend.box = element_blank(),
        legend.spacing.x = unit(1, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

dist_ecdf = ggplot(data = recovery_df %>% st_drop_geometry() %>% group_by(release_year, years_at_liberty, years_at_liberty_num) %>% summarise(recoveries = n()) %>% 
         group_by(years_at_liberty) %>% mutate(recoveries = recoveries / max(recoveries)), 
       aes(x = recoveries, color = years_at_liberty)) +
  stat_ecdf(lwd = 1.3, alpha = 0.85) +
  scale_color_viridis_d(option = 'magma') +
  labs(x = "Relative Number of Recoveries", y = 'Empirical Cumulative Distribution Function', col = "Time at Liberty (Years)", fill = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.85, 0.25), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank(), legend.box = element_blank(),
        legend.spacing.x = unit(1, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

dist_ecdf

### Distance Traveled ------------------------------------------------------
# get distance travelled
dist_tag_recovery_df = recovery_df %>% 
  mutate(dist_trav = distHaversine(cbind(HLNG, HLAT), cbind(RLNG, RLAT)))

dist_trav = ggplot(dist_tag_recovery_df, aes(x = dist_trav * 1e-3)) +
  geom_histogram(alpha = 0.5, col = 'black') +
  labs(x = "Distance Travelled (km)", y = "Count") +
  theme_bw(base_size = 15) +
  facet_wrap(~years_at_liberty, scales = 'free_y', ncol = 2) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

dist_trav

### Empirical Movement ------------------------------------------------------

## Bring in tag data
tag_recovery_df = readRDS(file = file.path("Data", "Tag_recoveryDF.RDS")) %>% st_drop_geometry()
tag_release_df = readRDS(file = file.path("Data", "Tag_releaseDF.RDS")) %>% st_drop_geometry()

##### Residual Munging --------------------------------------------------------
## separate inner-gulf releases CH and CL
area_labs = data.frame(actual_lab = c("WGOA","Inner-EGOA","CGOA","WGOA", "CGOA", "WY","EY/SE","BS","AI"),
                       obs_catch_lab = c("WGOA", "EGOA (outer & innner)","CGOA","Western Gulf of Alaska", "Central Gulf of Alaska", "West Yakutat", "East Yakutat / Southeast Alaska", "Bering Sea" , "Aleutian Islands" ))

tag_recovery_df$area_lab = area_labs$actual_lab[match(tag_recovery_df$NAME, area_labs$obs_catch_lab)]
tag_release_df$area_lab = area_labs$actual_lab[match(tag_release_df$NAME, area_labs$obs_catch_lab)]

areas_of_interest = c("BS","AI", "WGOA", "CGOA","WY", "EY/SE")
tag_recovery_df = tag_recovery_df %>% filter(area_lab %in% areas_of_interest, !is.na(area_lab))
tag_release_df = tag_release_df %>% filter(area_lab %in% areas_of_interest, !is.na(area_lab))

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
init_records = nrow(tag_recovery_df)
tag_recovery_df = tag_recovery_df %>% filter(TAG_TYPE.y %in% tag_types_allowed)
init_records = nrow(tag_release_df)
tag_release_df = tag_release_df %>% filter(TAG_TYPE %in% tag_types_allowed)

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

survey_age_data = readRDS(file = file.path("Data",  "Survey", "age_df_for_AF.RDS")) # created in DataCharacterisation/SurveyData.R

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

temp = tag_recovery_df %>% inner_join(tag_release_df, by = "recovery_id")
## when we merge area_lab.x is recovery region and area_lab.y is release region
temp_release = tag_release_df %>% full_join(tag_recovery_df, by = "recovery_id")
## when we merge area_lab.x is RELEASE region and area_lab.y is RECOVERY region
temp_release = temp_release %>% group_by(area_lab.y, area_lab.x) %>% summarise(releases = n(), recoveries = sum(!is.na(TAG_TYPE.x)))
temp_release$area_lab.y = factor(temp_release$area_lab.y, levels = rev(c("BS", "AI","WGOA","CGOA","WY", "EY/SE")), ordered = T)
temp_release$area_lab.x = factor(temp_release$area_lab.x, levels = rev(c("BS", "AI","WGOA","CGOA","WY", "EY/SE")), ordered = T)
sum_tag_release_df = temp_release %>% group_by(area_lab.x) %>% mutate(prop = recoveries / sum(recoveries), releases= sum(releases))
# drop NA regions
sum_tag_release_df = sum_tag_release_df %>% filter(!is.na(area_lab.y))

sum_tag_release_df$area_lab_release = paste0(sum_tag_release_df$area_lab.x, "\n(", sum_tag_release_df$releases,")")
unique(sum_tag_release_df$area_lab_release)
sum_tag_release_df$area_lab_release = factor(sum_tag_release_df$area_lab_release, levels = rev(c(    "BS\n(28126)" , "AI\n(18917)","WGOA\n(30875)",   "CGOA\n(81164)",  "WY\n(37969)",  "EY/SE\n(79373)")), ordered = T)


##### Plot --------------------------------------------------------------------
emp_move = ggplot(sum_tag_release_df, aes(x = area_lab.y, y = area_lab_release, fill = prop)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma', alpha = 0.5) +
  geom_text(aes(x = factor(area_lab.y), y = factor(area_lab_release), label = round(recoveries,2)), color = "black", size = 4) +
  labs(x = "Recovery Region", y = "Release Region (Number of releases)", fill = "Proportion of\nrecoveries") +
  theme_test(base_size = 15) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA)) 

emp_move

# Proportion of catch vs. samples -----------------------------------------
observer_length_df = read.csv(file.path("Data", "Observer", "norpac_length_report.csv"), skip = 6)
observer_age_df = read.csv(file.path("Data", "Observer", "norpac_age_report.csv"), skip = 6) 
observer_catch_df = read.csv(file.path("Data", "Observer", "norpac_catch_report.csv"), skip = 6)
longline_areas = st_read(file.path("Data", "managementareas", "Sablefish_Longline_Area","Sablefish_Longline_Area.shp"))
catch_post_90 = read.csv(file = file.path("Data", "Raw", "raw_catch_file.csv"))

### Residual Munging --------------------------------------------------------
# add catch data to this
observer_length_df_full = observer_length_df %>% left_join(observer_catch_df, by = "Haul.Join")
observer_age_df_full = observer_age_df %>% left_join(observer_catch_df, by = "Haul.Join")

observer_catch_df = observer_catch_df %>% mutate(fmp_gear =   
                                                  case_when(Gear.Description == "NON PELAGIC"  ~ "TRW",
                                                            Gear.Description == "PELAGIC"  ~ "TRW",
                                                            Gear.Description == "PAIR TRAWL"  ~ "TRW",
                                                            Gear.Description == "POT OR TRAP"  ~ "POT",
                                                            Gear.Description == "LONGLINER"  ~ "HAL"))

observer_catch_sf = st_as_sf(x = observer_catch_df,                         
                             coords = c("Lon.DD.End", "Lat.DD.End"),
                             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
observer_catch_nad83_df = st_transform(observer_catch_sf, crs = st_crs(longline_areas))
catch_nad83_df = st_transform(observer_catch_sf, crs = st_crs(longline_areas))
catch_nad83_df = catch_nad83_df %>% st_join(longline_areas)

obs_catch_labs = data.frame(actual_lab = c("WGOA", "CGOA", "WY","EY/SE","BS","AI"),
                            obs_catch_lab = c("Western Gulf of Alaska", "Central Gulf of Alaska", "West Yakutat", "East Yakutat / Southeast Alaska", "Bering Sea" , "Aleutian Islands" ))
catch_labs = data.frame(actual_lab = c("WGOA", "CGOA", "WY","EY/SE","EY/SE","BS","AI"),
                        catch_lab = c("WG", "CG", "WY", "EY", "SE" ,"BS" , "AI" ))

catch_nad83_df$area_lab = obs_catch_labs$actual_lab[match(catch_nad83_df$NAME, obs_catch_labs$obs_catch_lab)]
catch_post_90$area_lab = catch_labs$actual_lab[match(catch_post_90$fmp_subarea, catch_labs$catch_lab)]

catch_by_year_region = catch_post_90 %>% group_by(year, area_lab) %>% summarise(catch = sum(weight_posted)) %>% rename(Year = year)
# catch_by_year_region = catch_by_year_region %>% mutate(area_lab = case_when(area_lab == "EY"  ~ "EGOA",
#                                                                             area_lab == "WY"  ~ "EGOA",
#                                                                             TRUE ~ area_lab))

# summarise by year and area
reported_catch_by_year = catch_post_90 %>% group_by(area_lab, year) %>% summarise(catch_post_90 = sum(weight_posted))

obs_catch_by_year = as.data.frame(catch_nad83_df) %>% group_by(area_lab, Year) %>% summarise(catch_post_90 = sum( Extrapolated.Weight..kg.) / 1000)#, observed_events = length(Official.Total.Catch..mt.))
colnames(obs_catch_by_year) = colnames(reported_catch_by_year)
obs_catch_by_year$type = "Observed Catch"
reported_catch_by_year$type = "Reported Samples"
full_catch_by_year = rbind(reported_catch_by_year, obs_catch_by_year)
# calculate porportion of catch by method in a year, and see how well it was observed 
full_catch_by_year = full_catch_by_year %>% group_by(year, type) %>% mutate(percentage_catch = catch_post_90 / sum(catch_post_90) * 100)

### Plot --------------------------------------------------------------------
## summarise by area and year
shpe_manual = c("Reported Samples" = 1, "Observed Catch" = 3)
col_manual = c("Reported Samples" = "red", "Observed Catch" = "blue")

prop_fish_catch_effort_samp = ggplot(full_catch_by_year %>% filter(year>= 1990) %>% 
                                       mutate(area_lab = factor(area_lab, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE"))),
                                     aes(y = year, x = area_lab, col = type, shape = type)) +
  geom_point(aes(size = percentage_catch * if_else(type == "Observed Catch", 0.5, 1))) +
  scale_shape_manual(values=shpe_manual) +
  scale_size_area(max_size = 10, labels = c(10, 20), breaks = c(10, 20)) +
  scale_color_manual(values = col_manual) +
  scale_y_reverse() +
  labs(x = "Area", y = "Year", size  = "Percentage of Catch in a Year", col = "", shape = "", title = "") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_discrete(na.translate = FALSE) 

prop_fish_catch_effort_samp

# Abundance Indices and Catch ---------------------------------------------
full_catch_df = readRDS(file = here("Data", "5-Area", "Catch_by_year_area_gear.RDS"))
fixed_gear_with_imputation = readRDS(file = here("Data", "5-Area", "fixed_gear_with_imputations_S1.RDS"))
trawl_gear_with_imputation = readRDS(file = here("Data", "5-Area", "trawl_gear_with_imputations_S1.RDS"))
# srv_area = readRDS(here("Data", "Survey", "regional_abundance_estimates.RDS"))
srv_area = read.csv(here("Data", "Survey", "Area RPNs for Strata 3 to 7.csv"))

# convert to metric tonnes
full_catch_df$catch_mt = full_catch_df$catch / 1000
full_catch_df$type = "Non-Imputed Catch"

# Bind imputed catch
# fixed gear
fixed_gear_with_imputation = fixed_gear_with_imputation %>% 
  select(fmp_gear, year, area, imputed_catch) %>% 
  mutate(catch = imputed_catch * 1e3, type = 'Imputed Catch') %>% 
  rename(catch_mt = imputed_catch) %>% 
  filter(year < 1977)

# trawl gear
trawl_gear_with_imputation = trawl_gear_with_imputation %>% 
  select(fmp_gear, year, area, imputed_catch) %>% 
  mutate(catch = imputed_catch * 1e3, type = 'Imputed Catch') %>% 
  rename(catch_mt = imputed_catch) %>% 
  filter(year < 1977)

catch_area_gear = rbind(full_catch_df, fixed_gear_with_imputation, trawl_gear_with_imputation) %>% 
  mutate(fmp_gear = ifelse(fmp_gear == "HAL", "Fixed", "Trawl"))

# Catch
catch_area = ggplot(catch_area_gear %>% 
                      mutate(area = factor(area, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA'))), 
                    aes(x = year, y = catch_mt, lty = type)) +
  geom_line(lwd = 1) +
  facet_grid(fmp_gear ~ area) +
  labs(x = "Year", y = "Catch (mt)", lty = 'Type') +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "transparent", colour = NA)) 

# Survey 
srv_idx_area = ggplot(srv_area %>% 
                        mutate(Mgmt_Abrev = case_when(
                          NPFMC.Sablefish.Management.Area == "Aleutians" ~ "AI",
                          NPFMC.Sablefish.Management.Area == "Bering Sea" ~ "BS",
                          NPFMC.Sablefish.Management.Area == "Western Gulf of Alaska" ~ "WGOA",
                          NPFMC.Sablefish.Management.Area == "Central Gulf of Alaska" ~ "CGOA",
                          NPFMC.Sablefish.Management.Area == "West Yakutat" ~ "WY",
                          NPFMC.Sablefish.Management.Area == "East Yakutat/Southeast" ~ "EY/SE"
                        ),
                        Mgmt_Abrev = factor(Mgmt_Abrev, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE"))) %>% 
                        group_by(Country, Mgmt_Abrev, Year) %>% 
                        summarize(RPN = sum(RPN, na.rm = T)) %>% 
                        filter(Year <= 2021),
                      aes(x = Year, y = RPN, color = Country, fill = Country, shape = Country, lty = Country)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 4) +
  facet_wrap(~Mgmt_Abrev, nrow = 1) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = "Year", y = "Relative Population Numbers")

# Growth Differences ------------------------------------------------------
# Read in age data
age_dat_srv = read.csv(here('Data',  'Survey', 'age_view.csv'), check.names = FALSE)
age_dat_obs = read.csv(file.path("Data", "Observer", "norpac_age_report.csv"), skip = 6) 

# Clean up some age data
age_dat_srv = age_dat_srv %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`,
         Haul = `Station Number`,
         Mgmt_Area = `NPFMC Sablefish Mgmt Area`) %>% 
  filter(!Age %in% c(0, 1)) %>% 
  select(Sex_name, Length, Weight, Age, Year, Haul, Mgmt_Area) %>% 
  mutate(Sex_name = ifelse(Sex_name == "female", "Females", 'Males'),
         Mgmt_Area = case_when(
           Mgmt_Area == "Aleutians" ~ "AI",
           Mgmt_Area == "Bering Sea" ~ "BS",
           Mgmt_Area == "Western Gulf of Alaska" ~ "WGOA",
           Mgmt_Area == "Central Gulf of Alaska" ~ "CGOA",
           Mgmt_Area == "West Yakutat" ~ "WY",
           Mgmt_Area == "East Yakutat/Southeast" ~ "EY/SE"
           ), Type = "Survey")

# fishery
age_dat_obs = age_dat_obs %>% 
  filter(!is.na(Age), !is.na(`Length..cm.`), !is.na(`Weight..kg.`),
         Sex %in% c("F", "M"), !Age %in% c(0,1)) %>% 
  rename(Length = `Length..cm.`,
         Weight = `Weight..kg.`,
         Haul = Haul.Join,
         Mgmt_Area = FMP.Subarea,
         Sex_name = Sex) %>% 
  mutate(Gear_Abbrev = case_when(
    Gear.Description == "LONGLINER" ~ "Longliner",
    Gear.Description == "NON PELAGIC" ~ "Trawl",
    Gear.Description == "PELAGIC" ~ "Trawl",
    Gear.Description == "POT OR TRAP" ~ "Pot"
  ), Type = paste("Fishery", Gear_Abbrev)) %>% 
  select(Sex_name, Length, Weight, Age, Year, Haul, Mgmt_Area, Type) %>% 
  mutate(Sex_name = ifelse(Sex_name == "F", "Females", "Males"),
         Mgmt_Area = case_when(
           Mgmt_Area == "AI" ~ "AI",
           Mgmt_Area == "BS" ~ "BS",
           Mgmt_Area == "WG" ~ "WGOA",
           Mgmt_Area == "CG" ~ "CGOA",
           Mgmt_Area == "WY" ~ "WY",
           Mgmt_Area == "SE" ~ "EY/SE"
         ))

age_dat = rbind(age_dat_obs, age_dat_srv) %>% 
  mutate(Mgmt_Area = factor(Mgmt_Area, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'WY', "EY/SE")))

# Plot
grwth_area_gear = ggplot(age_dat %>% filter(!str_detect(Type, "Trawl")) %>% 
         mutate(), aes(x = Age, y = Weight, color = Mgmt_Area)) +
  geom_point(alpha = 0.05, size = 2.5) + 
  geom_smooth(method = 'gam', lwd = 1.3, formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_manual(values = c("black", "#56B4E9", "#009E73","#0072B2", "#D55E00", "#E69F00")) +
  facet_grid(Sex_name~Type, scales = 'free') +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        plot.background = element_rect(fill = "transparent", colour = NA))  +
  labs(x = "Age", y = "Weight (kg)", color = "AK Sablefish Mgmt Regions")

grwth_area_gear


# Manuscript Plots --------------------------------------------------------


### Data Comparison (1-Area vs. 5-Area) --------------------------------------------------------------------
# 1 Area model
ggsave(
  plot_obs(data = readRDS(here("Output", 'Final Models', "1-Area-1960", "Data.RDS")), 
           region_key = data.frame(area = c("1-Area"), TMB_ndx = c(0)), 
           survey_labels = c("Japanese LL Survey","Domestic LL survey")),
  filename = here("figs", "Manuscript_Plots", "1Area_Dat.png"),
  width = 17, height = 10
)

# 5 Area model
ggsave(
  plot_obs(readRDS(here("Output", 'Final Models', "5-Area-1960", "Data.RDS")), 
           region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4)), 
           survey_labels = c("Japanese LL Survey","Domestic LL survey")),
  filename = here("figs", "Manuscript_Plots", "5Area_Dat.png"),
  width = 17, height = 17
)


### Regression Tree Length Comps + Age Comps --------------------------------
ggsave(
  plot_grid(srv_reg_tree_map, fsh_reg_tree_map, srv_agg_ages, fish_agg_ages, 
            nrow = 2, rel_heights = c(0.6, 0.4), align = 'hv', axis = 'lr',
            labels = c('A', 'B', 'C', 'D'), label_size = 18, label_x = 0.025),
  filename = here("figs", "Manuscript_Plots", "Comp_Exploration.png"),
  width = 16, height = 10
)


### Recovery Map Plot -------------------------------------------------------
ggsave(
  recovery_map_plot,
  filename = here("figs", "Manuscript_Plots", "Recovery_Map.png"),
  width = 17, height = 10
)

### Cohort Age Map ----------------------------------------------------------
mean_cohort_dist = rbind(survey_cohort %>% mutate(obs_type = 'Survey'),
                         fish_cohort %>% mutate(obs_type = 'Fishery'))

# Survey
srv_cohort_map = ggplot() +
  geom_sf(data = west) +
  geom_sf(data = nmfs_areas, aes(lty = GEN_NAME), alpha = 0.5, lwd = 0.4) +
  geom_path(mean_cohort_dist %>% filter(Cohort %in% c(seq(1996, 2016, 2)), obs_type == 'Survey'),
          mapping = aes(x = mean_lon, y = mean_lat, color = Year, group = 1), lwd = 1) +
  geom_point(mean_cohort_dist %>% filter(Cohort %in% c(seq(1996, 2016, 2)), obs_type == 'Survey'),
             mapping = aes(x = mean_lon, y = mean_lat, color = Year, fill = Year, group = 1), colour="black",pch = 21, size = 5, alpha = 0.75) +
  coord_sf(ylim = c(48, 70.5), xlim = c(165, 227)) + # Restrict Map Area
  scale_fill_viridis_c(option = 'magma') +
  scale_color_viridis_c(option = 'magma') +
  facet_wrap(~Cohort, nrow = 2) +
  theme_bw(base_size = 11) +
  guides(lty = 'none') +
  theme(legend.position = 'right', legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12)) +
  labs(x = 'Longitude', y = 'Latitude')

# Fishery
fsh_cohort_map = ggplot() +
  geom_sf(data = west) +
  geom_sf(data = nmfs_areas, aes(lty = GEN_NAME), alpha = 0.5, lwd = 0.4) +
  geom_path(mean_cohort_dist %>% filter(Cohort %in% c(seq(1996, 2016, 2)), obs_type == 'Fishery'),
            mapping = aes(x = mean_lon, y = mean_lat, color = Year, group = 1), lwd = 1) +
  geom_point(mean_cohort_dist %>% filter(Cohort %in% c(seq(1996, 2016, 2)), obs_type == 'Fishery'),
             mapping = aes(x = mean_lon, y = mean_lat, color = Year, fill = Year, group = 1), colour="black",pch = 21, size = 5, alpha = 0.75) +
  coord_sf(ylim = c(48, 70.5), xlim = c(165, 227)) + # Restrict Map Area
  scale_fill_viridis_c(option = 'magma') +
  scale_color_viridis_c(option = 'magma') +
  facet_wrap(~Cohort, nrow = 2) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'none', legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12)) +
  labs(x = 'Longitude', y = 'Latitude')

# Get legend
legend = get_legend(srv_cohort_map)

# combine plots
cohort_map_plots = plot_grid(srv_cohort_map + theme(legend.position = 'none'), fsh_cohort_map, nrow = 2, align = 'hv', axis = 'lr',
                             labels = c('A', 'B'), label_size = 18, label_x = 0.03)
# save
ggsave(plot_grid(cohort_map_plots, legend, ncol = 2, rel_widths = c(1, 0.05)),
  filename = here("figs", "Manuscript_Plots", "Cohort_Comp_Exploration.png"),
  width = 16, height = 10
)


### Tag Recovery + Liberty + Distribution Plot ------------------------------
tag_recovery_liberty_plot = plot_grid(tag_recovery_ts, liberty_hist, nrow = 2,
                                      labels = c('A', 'B'), label_size = 15, label_x = 0.015)

ggsave(
  plot_grid(tag_recovery_liberty_plot, dist_ecdf, ncol = 2, rel_widths = c(0.45, 0.55),
            labels = c('', 'C'), label_size = 15, label_x = 0.015),
  filename = here("figs", "Manuscript_Plots", "Tag_Data_Exploration.png"),
  width = 17, height = 10
)


### Tag Distance + Empirical Movement ---------------------------------------

# distance traveled
ggsave(
  dist_trav,
  filename = here("figs", "Manuscript_Plots", "Tag_Dist_Trav.png"),
  width = 10, height = 13,
)

# enmpirical movement
ggsave(
  emp_move,
  filename = here("figs", "Manuscript_Plots", "Tag_Emp_Move.png"),
  width = 10, height = 7,
)


### Survey and Fishery effort by depth, longitude, latitude ----------------------------------------------------

# depth for survey
ggsave(
  srv_depth_len_plot,
  filename = here("figs", "Manuscript_Plots", "Srv_Len_Depth.png"),
  width = 16, height = 10,
)

# depth for fishery
ggsave(
  bot_depth_len_obs,
  filename = here("figs", "Manuscript_Plots", "Fsh_Len_Depth.png"),
  width = 16, height = 10,
)

# longitude for fishery
ggsave(
  lon_len_obs,
  filename = here("figs", "Manuscript_Plots", "Fsh_Len_long.png"),
  width = 16, height = 10,
)

# latitude for fishery
ggsave(
  lat_len_obs,
  filename = here("figs", "Manuscript_Plots", "Fsh_Len_lat.png"),
  width = 16, height = 10,
)


### Proportion of sampling effort vs. catch ---------------------------------

ggsave(
  prop_fish_catch_effort_samp,
  filename = here("figs", "Manuscript_Plots", "Fsh_Cat_Samp_Effort_prop.png"),
  width = 16, height = 10
)


### Catch and Survey Indicies -----------------------------------------------

# Catch
ggsave(
  catch_area,
  filename = here("figs", "Manuscript_Plots", "Fsh_Cat_Area.png"),
  width = 16, height = 10
)

# Survey Index
ggsave(
  srv_idx_area,
  filename = here("figs", "Manuscript_Plots", "Srv_Idx_Area.png"),
  width = 16, height = 10
)


### Growth Differences ------------------------------------------------------

ggsave(
  grwth_area_gear,
  filename = here("figs", "Manuscript_Plots", "Grwth_Diff_Area.png"),
  width = 16, height = 10
)



### Base Map ----------------------------------------------------------------

# Save map for making a conceptual model 
base_map = ggplot() +
  geom_sf(data = west, lwd = 0.2, color = 'black') + # World Map
  geom_sf(data = nmfs_areas, aes(lty = GEN_NAME), fill = NA, alpha = 0.35, lwd = 0.4) +
  
  # # Annotate with labels
  # annotate("text", label = "Alaska", x = 208, y = 65, size = 8) + # Alaska label
  # annotate("text", label = "Canada", x = 224.65, y = 63, size = 8) + # Canada label
  # annotate("text", label = "Russia", x = 172, y = 67, size = 8) + # Canada label
  
  coord_sf(ylim = c(45, 70.5), xlim = c(165, 235)) + # Restrict Map Area
  theme_bw(base_size = 15) +
  theme(legend.position = 'none', legend.box = "vertical", legend.background = element_blank(),
        legend.spacing.y = unit(0.075, "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 23),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = "Longitude", y = "Latitude")

ggsave(
  base_map,
  filename = here("figs", "Manuscript_Plots", "Base_Map.png"),
  width = 20, height = 15
)
