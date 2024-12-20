# Purpose: To plot spatial model final results for manuscript
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 10/13/24


# Set up ------------------------------------------------------------------
library(ggsflabel)
library(here)
library(tidyverse)
library(FishFreqTree)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(geosphere)
library(SpatialSablefishAssessment)
library(expm)
library(scico)

# Get model results
path = here("Output", "Final Models", "5-Area-1960")
mod_path = here(path, "5-Area-1960-03-FishBlock")
model_rep = readRDS(file.path(paste(mod_path, "/mle_report.RDS", sep = ""))) # model model report
model_5area_dat = readRDS(file.path(paste(mod_path, "/data.RDS", sep = ""))) # model model report

# Read in maps here
west = ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west = st_shift_longitude(west) # shift ongitude for plotting

# Read in stat areas
nmfs_areas = read_sf(dsn = here("data", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas = nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>% 
  mutate(         NAME = case_when(
    NAME == "Aleutian Islands" ~ "AI",
    NAME == "Bering Sea" ~ "BS",
    NAME == "Western Gulf of Alaska" ~ "WGOA",
    NAME == "Central Gulf of Alaska" ~ "CGOA",
    NAME == "West Yakutat" ~ "EGOA",
    NAME == "East Yakutat / Southeast Alaska" ~ "EGOA"
  ), NAME = factor(NAME, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>% 
  group_by(NAME) %>%
  summarise(geometry = st_union(geometry))

# Coerce longline areas
nmfs_areas = st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas = nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas = st_shift_longitude(nmfs_areas) # shift longitude for plotting

# get centroids of the geometry for plotting
centroids = nmfs_areas %>%
  group_by(NAME) %>%
  summarise(geometry = st_centroid(geometry)) %>%
  ungroup()

# Movement Dynamics -------------------------------------------------------
colors = unname(ggthemes::ggthemes_data[["colorblind"]][["value"]]) # get colors
### Residual Munging --------------------------------------------------------
move = model_rep$movement_matrix # get movement matrix

# Get Stationary Movement Dynamics
stationary_dist = data.frame()
time_seq = seq(1,25,1)
for(a in 1:dim(move)[4]) {
  for(t in 1:length(time_seq)) {
    stationary_dist_mat = move[,,1,a] %^% time_seq[t] # get stationary distribution
    colnames(stationary_dist_mat) = c("BS", "AI", "WGOA", "CGOA", "EGOA") # do some naming
    # munge into dataframe
    stationary_dist_tmp = data.frame(From = c("BS", "AI", "WGOA", "CGOA", "EGOA"), 
                                     stationary_dist_mat, 
                                     AgeBlock = paste("AgeBlock", a),
                                     Time = time_seq[t])
    stationary_dist = rbind(stationary_dist, stationary_dist_tmp)
  } # end t
}

# relevel
stationary_dist = stationary_dist %>% 
  pivot_longer(!c("From", "Time", "AgeBlock"), names_to = "To", values_to = "move") %>% 
  mutate(From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
         To = factor(To, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# rename variables
move_df = reshape2::melt(move) %>% 
  mutate(    Var1 = case_when(
    Var1 == 1 ~ "BS",
    Var1 == 2 ~ "AI",
    Var1 == 3 ~ "WGOA",
    Var1 == 4 ~ "CGOA",
    Var1 == 5 ~ "EGOA"
  ), Var1 = factor(Var1, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
  Var2 = case_when(
    Var2 == 1 ~ "BS",
    Var2 == 2 ~ "AI",
    Var2 == 3 ~ "WGOA",
    Var2 == 4 ~ "CGOA",
    Var2 == 5 ~ "EGOA"
  ), Var2 = factor(Var2, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>% 
  # filter(Var4 == 1) %>% 
  dplyr::select(Var1, Var2, Var4, value) %>% 
  rename(From = Var1, To = Var2, AgeBlk = Var4)

# Left join centroids to dataframe for plotting where from -> to
move_df_with_coords = move_df %>%
  left_join(centroids, by = c("From" = "NAME")) %>%
  rename(from_geometry = geometry) %>%
  left_join(centroids, by = c("To" = "NAME")) %>%
  rename(to_geometry = geometry) %>%
  mutate(
    from_x = st_coordinates(st_shift_longitude(from_geometry))[, 1], # get from centroids
    from_y = st_coordinates(from_geometry)[, 2],
    to_x = st_coordinates(st_shift_longitude(to_geometry))[, 1], # get to centroids
    to_y = st_coordinates(to_geometry)[, 2],
    AgeBlk = paste('AgeBlock', AgeBlk)
  ) 

### Movement Map ------------------------------------------------------------
map_movement_plot = ggplot() +
  geom_sf(data = nmfs_areas, alpha = 0.55) +
  geom_sf(data = west, lwd = 0.05, color = 'black', alpha = 1) + # World Map
  geom_curve(data = move_df_with_coords %>% filter(From != To),
             aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = To, size = value),
             curvature = 0.4, alpha = 0.45) + # movement arrows (from != to)
  geom_text(data = move_df_with_coords %>% filter(From != To),
             aes(x = to_x, y = to_y, label = round(value, 2), color = To), alpha = 1, size = 8.5, nudge_y = 0.07) + # from != to
  scale_color_manual(values = colors[-c(1,5)]) +
  geom_text(data = move_df_with_coords %>% filter(From == To), size = 8.5, # from == to 
             aes(x = from_x, y = from_y, label = round(value, 2)), alpha = 1, color = 'black', nudge_y = 0.07) +
  facet_grid(AgeBlk~From) +
  coord_sf(ylim = c(45.2, 70.5), xlim = c(165, 230)) + # Restrict Map Area
  theme_bw(base_size = 24) +
  labs(x = "Longitude", y = "Latitude", size = 'P(movement)',
       color = "To", label = "P(residence)") +
  theme(legend.position = 'top',
        legend.box = 'vertical',
        legend.background = element_blank(),
        legend.spacing.y = unit(-1, 'cm'),  # Adjust spacing between items
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

### Stationary Movement -------------------------------------------------
stationary_movement_plot = ggplot(stationary_dist, aes(x = To, y = move, color = Time, group = Time)) +
  geom_line(lwd = 1.3, alpha = 0.3) +
  geom_point(size = 6) + 
  facet_grid(AgeBlock~From) +
  scale_color_viridis_c(option = 'magma') +
  theme_bw(base_size = 24) +
  labs(x = "Region (To)", y = "Stationary Movement Probability") +
  theme(legend.position = c(0.87,0.2), legend.key.height=unit(0.5,"cm"), legend.key.width=unit(0.95,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1), 
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_blank(), legend.direction = 'horizontal')

ggsave(
  plot_grid(map_movement_plot, stationary_movement_plot, ncol = 1,
            labels = c("A", 'B'), label_size = 28, hjust = -0.5, rel_heights = c(0.6, 0.4)), 
  filename = here("figs", "Manuscript_Plots", "Movement_Dynamics.png"),
  width = 25, height = 23,
) 

# Time Series Comparison --------------------------------------------------
# Get spatial ssb and recruitment
spatial_ssb = reshape2::melt(model_rep$SSB_yr) %>% mutate(type = 'Sp Bio (kt)')
spatial_dep = reshape2::melt(apply(model_rep$SSB_yr, 1, function(x) x / model_rep$Bzero)) %>% mutate(type = 'Depletion') %>% rename(Var1 = Var2, Var2 = Var1)
spatial_rec = reshape2::melt(model_rep$recruitment_yr) %>% mutate(type = 'Recrmt (millions)')
spatial_mean_age = reshape2::melt(model_rep$natage_f) %>% mutate(type = 'Mean Female Age')
  
# Read in single-area models
sgl_area_rep = readRDS(file.path(paste(here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final"), "/mle_report.RDS", sep = ""))) # model model report
sgl_area_dat = readRDS(here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final", "data.RDS")) # single area data

# get single area ssb and recruitment
sgl_ssb = reshape2::melt(sgl_area_rep$SSB_yr) %>% mutate(type = 'Sp Bio (kt)')
sgl_dep = reshape2::melt(sgl_area_rep$SSB_yr / sgl_area_rep$Bzero) %>% mutate(type = 'Depletion')
sgl_rec = reshape2::melt(sgl_area_rep$recruitment_yr) %>% mutate(type = 'Recrmt (millions)')

### Residual Munging --------------------------------------------------------
# Munging for spatial mean female age
spatial_mean_age_sum = spatial_mean_age %>% 
  group_by(Var3, Var2, type) %>% # group by region and year
  mutate(prop = value / sum(value)) %>% 
  summarize(value = sum(prop * 2:31)) %>% 
  filter(Var3 != 63) %>% 
  rename(Var1 = Var3)

# get spatial time series
spatial_ts = rbind(spatial_ssb, spatial_dep, spatial_rec, spatial_mean_age_sum) %>% 
  rename(Year = Var1, region = Var2) %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), 
    region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    Year = Year + 1959) %>% 
  mutate(model = '5-Area')

# get total region wide time series
region_wide_ts = spatial_ts %>% 
  filter(!str_detect(type, 'Mean')) %>% # remove female age - do separately given slightly complex
  group_by(Year, type) %>% 
  summarize(value = sum(value)) %>% 
  mutate(model = '5-Area') 

# Replace depletion with alternative calc (not grouping and then summing)
region_wide_ts[region_wide_ts$type == 'Depletion',]$value = rowSums(model_rep$SSB_yr) / sum(model_rep$Bzero)

# get region wide mean age
region_wide_mean_age = spatial_mean_age %>% 
  mutate(Var1 = Var1 + 1) %>% 
  group_by(Var3, type) %>% # group by region and year
  mutate(prop = value / sum(value)) %>% 
  summarize(value = sum(prop * Var1)) %>% 
  filter(Var3 != 63) %>% 
  rename(Year = Var3) %>% 
  mutate(model = '5-Area',
         Year = Year + 1959)

# munge single area time series
sgl_ts = rbind(sgl_ssb, sgl_dep, sgl_rec) %>% 
  select(Var1, value, type) %>% 
  rename(Year = Var1) %>% 
  mutate(Year = Year + 1959, 
         model = '1-Area')

# single area get mean age
sgl_mean_age_f = reshape2::melt(sgl_area_rep$natage_f) %>% 
  rename(age = Var1, region = Var2, Year = Var3) %>% 
  group_by(Year) %>% 
  mutate(prop = value / sum(value)) %>% 
  summarize(value = sum(prop * 2:31)) %>% 
  mutate(Year = Year + 1959,
         type = 'Mean Female Age',
         model = '1-Area') %>% 
  filter(Year != 2022) 


#### Harvest Rate Munging ----------------------------------------------------
# single area get harvest rate
# Get weight at age to calcualte harvest rate
waa_f = array(sgl_area_dat$female_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), sgl_area_dat$n_regions, length(sgl_area_dat$years) + 1))
waa_m = array(sgl_area_dat$male_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), sgl_area_dat$n_regions, length(sgl_area_dat$years) + 1))

# Get array
# females
sel_fixed_f = aperm(replicate(1, cbind(replicate(length(1960:2015), sgl_area_rep$sel_fixed_f[,1]),
                                       replicate(length(2016:2022), sgl_area_rep$sel_fixed_f[,2]))), c(1,3,2))

# males
sel_fixed_m = aperm(replicate(1, cbind(replicate(length(1960:2015), sgl_area_rep$sel_fixed_m[,1]),
                                       replicate(length(2016:2022), sgl_area_rep$sel_fixed_m[,2]))), c(1,3,2))

# 1 Area harvest rate
# Get exploitable biomass for fixed gear at age
expl_fixed_age1 = sel_fixed_f * sgl_area_rep$natage_f * waa_f +
  sel_fixed_m * sgl_area_rep$natage_m * waa_m 
pred_fixed_cat1 = sgl_area_rep$fixed_fishery_catch # get predicted fixed gear catch
# Get exploitable biomass for trawl gear at age
expl_trwl_age1 = as.vector(sgl_area_rep$sel_trwl_f) * sgl_area_rep$natage_f * waa_f +
  as.vector(sgl_area_rep$sel_trwl_m) * sgl_area_rep$natage_m * waa_m 
pred_trawl_cat1 = sgl_area_rep$trwl_fishery_catch # get predicted fixed gear catch

# Sum across ages
expl_fixed1 = apply(expl_fixed_age1, c(2,3), sum)
expl_trwl1 = apply(expl_trwl_age1, c(2,3), sum)

# Single area harvest rate
sgl_hr_fixed = data.frame(Year = 1960:2021, type = 'Fixed-Gear Hrv Rate', value = t(pred_fixed_cat1 / expl_fixed1[-63]), model = '1-Area')
sgl_hr_trwl = data.frame(Year = 1960:2021, type = 'Trawl Gear Hrv Rate', value = t(pred_trawl_cat1 / expl_trwl1[-63]), model = '1-Area')


# Region-wide combined 5-area harvest rate
waa_f = array(model_5area_dat$female_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), model_5area_dat$n_regions, length(model_5area_dat$years) + 1))
waa_m = array(model_5area_dat$male_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), model_5area_dat$n_regions, length(model_5area_dat$years) + 1))

# females
sel_fixed_f = aperm(replicate(5, cbind(replicate(length(1960:2015), model_rep$sel_fixed_f[,1]),
                                       replicate(length(2016:2022), model_rep$sel_fixed_f[,2]))), c(1,3,2))

# males
sel_fixed_m = aperm(replicate(5, cbind(replicate(length(1960:2015), model_rep$sel_fixed_m[,1]),
                                       replicate(length(2016:2022), model_rep$sel_fixed_m[,2]))), c(1,3,2))

expl_fixed_age5 = sel_fixed_f * model_rep$natage_f * waa_f +
  sel_fixed_m * model_rep$natage_m * waa_m 
pred_fixed_cat5 = model_rep$fixed_fishery_catch # get predicted fixed gear catch
# Get exploitable biomass for trawl gear at age
expl_trwl_age5 = as.vector(model_rep$sel_trwl_f) * model_rep$natage_f * waa_f +
  as.vector(model_rep$sel_trwl_m) * model_rep$natage_m * waa_m 
pred_trawl_cat5 = model_rep$trwl_fishery_catch # get predicted fixed gear catch
# Sum across ages
expl_fixed5 = apply(expl_fixed_age5, c(2,3), sum)
expl_trwl5 = apply(expl_trwl_age5, c(2,3), sum)

# F area fixed-gear harvest rate
region_wide_hr_fixed = data.frame(Year = 1960:2021, type = 'Fixed-Gear Hrv Rate', value = colSums(pred_fixed_cat5) / colSums(expl_fixed5[,-63]), model = '5-Area')
region_wide_hr_trwl = data.frame(Year = 1960:2021, type = 'Trawl Gear Hrv Rate', value = colSums(pred_trawl_cat5) / colSums(expl_trwl5[,-63]), model = '5-Area')

# get spatially explicit harvest rate for fixed gear
spatial_hr_fixed = reshape2::melt(pred_fixed_cat5 / expl_fixed5[,-63]) %>% 
  rename(Year = Var2, region = Var1) %>% 
  mutate(type = 'Fixed-Gear Hrv Rate', model = '5-Area', Year = Year + 1959,
         region = case_when(
           region == 1 ~ "BS",
           region == 2 ~ "AI",
           region == 3 ~ "WGOA",
           region == 4 ~ "CGOA",
           region == 5 ~ "EGOA"
         ),  
         region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# get soatuakky explicit harvest rate for trawl gear
spatial_trwl_fixed = reshape2::melt(pred_trawl_cat5 / expl_trwl5[,-63]) %>% 
  rename(Year = Var2, region = Var1) %>% 
  mutate(type = 'Trawl Gear Hrv Rate', model = '5-Area', Year = Year + 1959,
         region = case_when(
           region == 1 ~ "BS",
           region == 2 ~ "AI",
           region == 3 ~ "WGOA",
           region == 4 ~ "CGOA",
           region == 5 ~ "EGOA"
         ),  
         region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# bind together
spatial_ts = rbind(spatial_ts, spatial_hr_fixed, spatial_trwl_fixed) %>% 
  mutate(type = factor(type, levels = c("Fixed-Gear Hrv Rate", "Trawl Gear Hrv Rate", "Mean Female Age",
                                        "Recrmt (millions)", "Sp Bio (kt)", "Depletion")))

# bind together
region_ts = rbind(region_wide_ts, sgl_ts, sgl_mean_age_f, 
                  region_wide_mean_age, sgl_hr_fixed, sgl_hr_trwl, region_wide_hr_fixed, region_wide_hr_trwl) %>% 
  mutate(type = factor(type, levels = c("Fixed-Gear Hrv Rate", "Trawl Gear Hrv Rate", "Mean Female Age",
                                        "Recrmt (millions)", "Sp Bio (kt)", "Depletion")))


### Time-Series Plots -------------------------------------------------------
# spatial time series plot
spatial_ts_plot = ggplot(spatial_ts, aes(x = Year, y = value, color = factor(region))) +
  geom_line(lwd = 0.85, alpha = 0.8) + 
  scale_color_manual(values = colors[-c(1,5)]) +
  facet_grid(type~region, scales = 'free') +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  labs(x = 'Year', y = 'Value') +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = "transparent", colour = NA))

region_wide_plot = ggplot(region_ts, aes(x = Year, y = value, lty = model)) +
  geom_line(lwd = 0.85, alpha = 0.8) + 
  facet_grid(type~"Region-Wide", scales = 'free') +
  theme_bw(base_size = 15) +
  coord_cartesian(ylim = c(0,NA)) +
  labs(x = 'Year', y = 'Value', lty = 'Model') +
  theme(legend.position = c(0.8, 0.78),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  plot_grid(spatial_ts_plot, region_wide_plot, rel_widths = c(0.75, 0.25), labels = c('A', 'B'),
            label_size = 20, hjust = -0.5),
  filename = here("figs", "Manuscript_Plots", "TimeSeries_Comparison.png"),
  width = 15, height = 11.5
) 


### Apportionment Time-Series Munging -----------------------------------------------

# Read in survey proprotions (from 2023 assessment files via Dan)
dom_srv_prop = readxl::read_xlsx(here("Data", "lls_RPW_area_prop.xlsx")) %>% 
  filter(five_yr_avg != "NA") %>% 
  mutate(area = case_when(
    area == "Aleutians" ~ "AI",
    area == "Bering Sea" ~ "BS",
    area == "Western Gulf of Alaska" ~ "WGOA",
    area == "Central Gulf of Alaska" ~ "CGOA",
    area == "West Yakutat" ~ "EGOA",
    area == "East Yakutat/Southeast" ~ "EGOA"
  ), 
  area = factor(area, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>% 
  rename(region = area)

# quick munging for survey proportions
dom_srv_prop = dom_srv_prop %>% 
  group_by(year, region) %>% 
  summarize(prop = sum(as.numeric(prop)),
            five_yr_avg = sum(as.numeric(five_yr_avg)))


# single area 
# Get weight at age to survey exploitable biomass
waa_f = array(sgl_area_dat$female_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), sgl_area_dat$n_regions, length(sgl_area_dat$years) + 1))
waa_m = array(sgl_area_dat$male_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), sgl_area_dat$n_regions, length(sgl_area_dat$years) + 1))

# filter to years of interst
waa_f = waa_f[,,1994:2021 - 1959]
waa_m = waa_m[,,1994:2021 - 1959]

# Get selex
sel_srv_f = cbind(replicate(length(1994:2021), sgl_area_rep$sel_srv_f[,1,2]))
sel_srv_m = cbind(replicate(length(1994:2021), sgl_area_rep$sel_srv_m[,1,2]))

# get NAA
naa_f = sgl_area_rep$natage_f[,1,1994:2021 - 1959]
naa_m = sgl_area_rep$natage_m[,1,1994:2021 - 1959]

# 1 Area survey exploitable biomass
expl_srv_age1 = sel_srv_f * waa_f * naa_f +
                sel_srv_m * waa_m * naa_m

# Exploitable survey biomass
expl_srv_biom1 <- unique(sgl_area_rep$srv_q[,,2]) * colSums(expl_srv_age1) 
expl_srv_ts1 <- data.frame(Year = 1994:2021, value = expl_srv_biom1)

# five area 
# Get weight at age to survey exploitable biomass
waa_f = array(sgl_area_dat$female_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), 5, length(sgl_area_dat$years) + 1))
waa_m = array(sgl_area_dat$male_mean_weight_by_age, dim = c(length(sgl_area_dat$ages), 5, length(sgl_area_dat$years) + 1))

# filter to years of interst
waa_f = waa_f[,,1994:2021 - 1959]
waa_m = waa_m[,,1994:2021 - 1959]

# Get selex
sel_srv_f = aperm(replicate(5, cbind(replicate(length(1994:2021), sgl_area_rep$sel_srv_f[,1,2]))), c(1,3,2))
sel_srv_m = aperm(replicate(5, cbind(replicate(length(1994:2021), sgl_area_rep$sel_srv_m[,1,2]))), c(1,3,2))

# get NAA
naa_f = model_rep$natage_f[,,1994:2021 - 1959]
naa_m = model_rep$natage_f[,,1994:2021 - 1959]

# 1 Area survey exploitable biomass
expl_srv_age5 = sel_srv_f * waa_f * naa_f +
                sel_srv_m * waa_m * naa_m

# Exploitable survey biomass
expl_srv_biom5 <- unique(model_rep$srv_q[,,2]) * apply(expl_srv_age5, c(2,3), sum)
expl_srv_ts5 <- reshape2::melt(expl_srv_biom5) %>% rename(region = Var1, Year = Var2) %>% 
  mutate(Year = Year + 1993,
         region = case_when(
           region == 1 ~ "BS",
           region == 2 ~ "AI",
           region == 3 ~ "WGOA",
           region == 4 ~ "CGOA",
           region == 5 ~ "EGOA"
         ),  
         region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))
  
# get ssb derived from survey biomass proprotions for a single-area model
sgl_ssb_srv_prop_df = dom_srv_prop %>% 
  left_join(expl_srv_ts1, by = c("year" = 'Year')) %>% 
  mutate(prop_value = prop * value, five_yr_value = five_yr_avg * value) %>% 
  left_join(expl_srv_ts5 %>% # left join spatial model values for comparison
              rename(five_area_value = value) %>% 
              select(Year, region, five_area_value), by = c("year" = "Year", "region")) %>% 
  drop_na() %>% group_by(year) %>% 
  mutate(five_area_prop = five_area_value / sum(five_area_value)) # just get the raw proportions

plot(expl_srv_biom1)
lines(colSums(expl_srv_biom5))

### Apportionment Time-Series Plot --------------------------------------------------------------------
apport_ts_plot = ggplot() +
  geom_line(sgl_ssb_srv_prop_df, 
            mapping = aes(x = year, y = five_yr_value, color = "1A Srv Idx Biom * Srv Prop (5 Yr Avg)"), lwd = 1.5, alpha = 0.85) + # 5-year average proportion
  geom_line(sgl_ssb_srv_prop_df, 
            mapping = aes(x = year, y = five_area_value, color = "5A Srv Idx Biom"), lwd = 1.5, alpha = 0.85) + # 5 area model
  facet_wrap(~region) +
  scale_color_manual(values = c("#21918c", "#3b528b")) +
  facet_grid(~region, scales = 'free') +
  theme_bw(base_size = 20) +
  coord_cartesian(ylim = c(0,NA)) +
  labs(x = 'Year', y = 'Srv Idx Biom',
       color = "", lty = "") +
  theme(legend.position = c(0.1, 0.93), plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_blank(), legend.text = element_text(size = 16))

# plot by proportions
# apport_ts_prop_plot = ggplot() +
#   geom_line(sgl_ssb_srv_prop_df, 
#             mapping = aes(x = year, y = five_yr_avg, color = "Srv Prop (5 Yr Avg)"), lwd = 1.5, alpha = 0.85) + # 5-year average proportion
#   geom_line(sgl_ssb_srv_prop_df, 
#             mapping = aes(x = year, y = five_area_prop, color = "5A Srv Idx Biom Prop"), lwd = 1.5, alpha = 0.85) + # 5 area model
#   facet_wrap(~region) +
#   scale_color_manual(values = c("#21918c", "#3b528b"),
#     breaks = c("Srv Prop (5 Yr Avg)", "5A Srv Idx Biom Prop")
#   ) + 
#   facet_grid(~region, scales = 'free') +
#   theme_bw(base_size = 20) +
#   coord_cartesian(ylim = c(0,NA)) +
#   labs(x = 'Year', y = 'Biom Prop',
#        color = "", lty = "") +
#   theme(legend.position = c(0.07, 0.93), plot.background = element_rect(fill = "transparent", colour = NA),
#         legend.background = element_blank(), legend.text = element_text(size = 16))

# regress the two together by biomass
regression_ts_plot = ggplot(sgl_ssb_srv_prop_df, aes(x = five_area_value, y = five_yr_value)) +
  geom_point() + 
  ggpmisc::stat_poly_line() +
  ggpmisc::stat_poly_eq(size = 6, label.y = 0.065, label.x = 0.97) +
  facet_wrap(~region) +
  scale_color_manual(values = c("#21918c", "#3b528b")) +
  facet_grid(~region, scales = 'free') +
  theme_bw(base_size = 20) +
  coord_cartesian(ylim = c(0,NA)) +
  labs(y = '1A Srv Idx Biom * Srv Prop (5 Yr Avg)', x = '5A Srv Idx Biom') +
  theme(legend.position = 'top', plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_blank())

# # regress the two together by proportions
# regression_prop_plot = ggplot(sgl_ssb_srv_prop_df, aes(x = five_area_prop, y = five_yr_avg)) +
#   geom_point() + 
#   ggpmisc::stat_poly_line() +
#   ggpmisc::stat_poly_eq(size = 6, label.y = 0.065, label.x = 0.97) +
#   facet_wrap(~region) +
#   scale_color_manual(values = c("#21918c", "#3b528b")) +
#   facet_grid(~region, scales = 'free') +
#   theme_bw(base_size = 20) +
#   coord_cartesian(ylim = c(0,NA)) +
#   labs(y = 'Srv Prop (5 Year Avg)', x = '5A Srv Idx Biom Prop') +
#   theme(legend.position = 'top', plot.background = element_rect(fill = "transparent", colour = NA),
#         legend.background = element_blank())


ggsave(
  plot_grid(apport_ts_plot, regression_ts_plot, ncol = 1,
            labels = c('A', 'B'), align = 'hv',
            label_size = 25, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "Apportionment_Comparison.png"),
  width = 23, height = 13
) 



# ISS Sensitivity ---------------------------------------------------------
# Look at ISS Sensitivity
out_5area_path = here("Output", "Final Models", "5-Area-1960")
files = list.files(out_5area_path)
files = files[str_detect(files, "FishBlock")]

ts_all <- data.frame()
for(i in 1:length(files)) {
  rep = readRDS(here(out_5area_path, files[i], 'mle_report.RDS')) # read in report
  mod = readRDS(here(out_5area_path, files[i], 'sd_report.RDS')) # read in report
  ssb_tmp <- reshape2::melt(rep$SSB_yr) %>% mutate(model = files[i], type = 'Sp Bio (kt)', conv = mod$pdHess) # get ssb
  rec_tmp <- reshape2::melt(rep$recruitment_yr) %>% mutate(model = files[i], type = 'Recrmt (millions)', conv = mod$pdHess) # get rec
  ts_all <- rbind(ts_all, ssb_tmp, rec_tmp)
} # end i

# do some more residual munging
ts_all <- ts_all %>% mutate(iss = str_remove(model, "5-Area-1960-03-FishBlock")) %>% 
  mutate(iss = as.numeric(ifelse(iss == "", 40, str_remove(iss, "_")))) %>% 
  mutate(Var2 = case_when(
    Var2 == 1 ~ "BS",
    Var2 == 2 ~ "AI",
    Var2 == 3 ~ "WGOA",
    Var2 == 4 ~ "CGOA",
    Var2 == 5 ~ "EGOA"
         ),  
    Var2 = factor(Var2, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# Plot!
ggplot(ts_all %>% filter(conv == TRUE), aes(x = Var1 + 1959, y = value, color = factor(iss))) +
  geom_line(lwd = 1.2) +
  scale_color_viridis_d(option = "cividis") +
  facet_grid(type~Var2, scales = "free_y") +
  labs(x = "Year", y = "Value", color = "Input Sample Size") +
  theme_bw(base_size = 15) 
