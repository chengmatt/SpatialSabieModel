#'
#' SurveyData.R
#' An exploratory analysis on the Longline survey
#'

# source("(00) Init.R")
library(tidyverse)
library(sf)
library(rgdal)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(RColorBrewer)
theme_set(theme_bw())
library(FishFreqTree) ## devtools::install_github('HaikunXu/RegressionTree',ref='main')
library(here)
## get some land objects for spatial plots
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## observed data projection decimal degrees, -180 -180 
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## read in management shape files
large_areas <- st_read(file.path("Data", "managementareas", "EGOA_CGOA_WGOA_AI","EGOA_CGOA_WGOA_AI.shp"))
longline_areas <- st_read(file.path("Data", "managementareas", "Sablefish_Longline_Area","Sablefish_Longline_Area.shp"))
st_crs(large_areas) ## look at projection coordinates
st_crs(longline_areas) ## look at projection coordinates
large_areas$Id = c("EGOA", "CGOA", "WGOA")
# choose a point on the surface of each geometry
large_area_points <- sf::st_point_on_surface(large_areas)
large_area_coords <- as.data.frame(sf::st_coordinates(large_area_points))
large_area_coords$NAME <- large_areas$Id
large_area_coords[2,2] = large_area_coords[2,2] - 1 ## drop the CGOA label slightly
large_area_coords[3,2] = large_area_coords[3,2] + 1 ## drop the CGOA label slightly
large_area_coords[3,1] = large_area_coords[3,1] - 2 ## drop the CGOA label slightly
LL_xlim = range(st_coordinates(longline_areas)[,"X"])
LL_ylim = range(st_coordinates(longline_areas)[,"Y"])
## LL areas
LL_area = st_area(longline_areas) * 1e-6 #km^2

## a data frame to match areas with two letter accronym
obs_catch_labs = data.frame(actual_lab = c("WGOA", "CGOA", "EGOA","EGOA","BS","AI"),
                            obs_catch_lab = c("Western Gulf of Alaska", "Central Gulf of Alaska", "West Yakutat", "East Yakutat / Southeast Alaska", "Bering Sea" , "Aleutian Islands" ))
catch_labs = data.frame(actual_lab = c("WGOA", "CGOA", "EGOA","EGOA","EGOA","BS","AI"),
                        catch_lab = c("WG", "CG", "WY", "EY", "SE" ,"BS" , "AI" ))


################
## Read in data
###############
survey_length_df = read.csv(file.path("Data", "Survey", "length_summary_view.csv")) ## this is observer data?
survey_age_df = read.csv(file.path("Data", "Survey", "age_view.csv")) ## this is observer data?
haul_df = read.csv(file.path("Data", "Survey", "Haul.csv"))
catch_df = read.csv(file.path("Data", "Survey", "catch_summary_view.csv"))
# - get derived indicies that are used in the assessment these are numbers
index_by_area_df = read.csv(file.path("Data", "Survey", "RPN_index_Strat_3_7.csv"))
index_df = read.csv(file.path("Data", "Survey", "Alaska_full_area.csv"))

## create an id to link among catch, age and length data sets
## id = Year-station.number-vesselnumber
catch_df$event_id = paste0(catch_df$Year, "-",catch_df$Station.Number, "-" ,catch_df$Vessel.Number)
survey_age_df$event_id = paste0(survey_age_df$Year, "-",survey_age_df$Station.Number, "-" ,survey_age_df$Vessel.Number)
survey_length_df$event_id = paste0(survey_length_df$Year, "-",survey_length_df$Station.Number, "-" ,survey_length_df$Vessel.Number)


## look at a single station to get comfortable with the data
test_station = catch_df %>% filter(Year == 1999, Station.Number == 80)

## is there a relationship between depth and catch
ggplot(catch_df, aes(x = Interpolated.Depth..m., y = Catch)) +
  geom_point(size = 0.3) +
  geom_smooth()

## aggregate data to station level i.e., over all skate/Hachi
station_catch = catch_df %>% group_by(event_id) %>% summarise(Year = unique(Year), Country = unique(Country), Vessel.Number = unique(Vessel.Number),    Vessel.Name = unique(Vessel.Name),
                                                              Station.Number = unique(Station.Number), NPFMC.Sablefish.Management.Area = unique(NPFMC.Sablefish.Management.Area), 
                                                              Habitat.Type = unique(Habitat.Type), Station.Type = unique(Station.Type), Geographic.Area.Name = unique(Geographic.Area.Name),
                                                              NMFS.Management.Area= unique(NMFS.Management.Area), Haul.Date = unique(Haul.Date), Haul = mean(Haul), Start.Latitude..DD. = mean(Start.Latitude..DD.),
                                                              Start.Longitude..DD. = mean(Start.Longitude..DD.), End.Latitude..DD. = mean(End.Latitude..DD.), End.Longitude..DD. = mean(End.Longitude..DD.), catch = sum(Catch),
                                                              mean_depth = mean(Interpolated.Depth..m.), start_time = unique(Haul.Start.Time..HHMM.), end_time = unique(Haul.End.Time..HHMM.), soad_time = unique(Soak.Time..min.),
                                                              gear_temp = unique(Gear.Temp..C.), surface_temp = unique(Surface.Temp..C.))


length(unique(station_catch$event_id))
length(unique(survey_age_df$event_id))
length(unique(survey_length_df$event_id))
## check all ids are 
table(unique(survey_length_df$event_id) %in% unique(station_catch$event_id))
table(unique(survey_age_df$event_id) %in% unique(station_catch$event_id))


## join catch data onto length data Cruise_Station_ID
length_df_full = survey_length_df %>% left_join(station_catch, by = "event_id")
age_df_full = survey_age_df %>% left_join(station_catch, by = "event_id")
## check dims after this
dim(survey_length_df)
dim(length_df_full)


################
## Characterise catch data first
###############
catch_sf = st_as_sf(x = station_catch,                         
                    coords = c("Start.Longitude..DD.", "Start.Latitude..DD."),
                    crs = projcrs)
## transform to more convenient projection "North American Datum 1983" 
## this removes the -180 problem
catch_nad83_df = st_transform(catch_sf, crs = st_crs(longline_areas))
catch_nad83_df = catch_nad83_df %>% st_join(longline_areas)
catch_nad83_df$area_lab = obs_catch_labs$actual_lab[match(catch_nad83_df$NAME, obs_catch_labs$obs_catch_lab)]

xlim = range(st_coordinates(catch_sf)[,"X"])
ylim = range(st_coordinates(catch_sf)[,"Y"])
if(FALSE) {
  ## plot catch data over all years
  all_station_plot = ggplot() +
    geom_sf(data = longline_areas, aes(col = "Longline"), fill = NA, linetype = "dashed", linewidth = 1) +
    geom_sf(data = world, fill = "darkgreen") +
    geom_sf(data = catch_sf, fill = "black") +
    coord_sf(crs  = st_crs(longline_areas), xlim = LL_xlim, ylim = LL_ylim) +
    labs(x = "", y = "", col = "") +
    theme(axis.text = element_blank(),
          legend.position = "bottom") +
    guides(col = "none")
  ggsave(plot = all_station_plot, filename = file.path(DIR$survey_figs, "survey_stations_all_years.png"), width = 8, height =4)
}

## rasterise catch over the years
## 50km by 50 km grid
raster_grd = raster(crs = crs(longline_areas), vals = NA, resolution = c(50000, 50000), ext = extent(c(LL_xlim, LL_ylim)))

years = unique(catch_sf$Year)
n_y = length(years)
z_breaks = c(seq(0,7, by = 1), 20)
z_labs = c(seq(0,6, by = 1),"7+", "")
display.brewer.pal(9, "RdYlBu")
z_pallete = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(length(z_breaks) - 1))
if(FALSE) {
  for(y in 1:n_y) {
    this_df = catch_nad83_df %>% filter(Year == years[y])
    this_rster  = rasterize(x = this_df, y = raster_grd, field = "catch", fun = sum, na.rm = T)
    this_rster_df = as.data.frame(this_rster, xy = TRUE)
    
    this_yr_catch = ggplot() +
      geom_sf(data = longline_areas, col = "gray60", fill = NA, linetype = "dashed", linewidth = 1) +
      geom_sf(data = world, fill = "darkgreen") +
      geom_raster(data = this_rster_df, aes(x = x, y = y, fill = layer / 1000)) +
      scale_fill_gradientn(na.value = NA, colours = z_pallete, breaks = z_breaks, limits = c(0,8), labels = z_labs) +
      coord_sf(crs  = st_crs(longline_areas), xlim = LL_xlim, ylim = LL_ylim) +
      labs(x = "", y = "", col = "", fill = "Numbers (000's)") +
      theme(axis.text = element_blank(),
            legend.position = "bottom") +
      guides(col = "none") +
      ggtitle(years[y])
    ggsave(plot = this_yr_catch, filename = file.path(DIR$survey_figs, "RasterCatch", paste0("catch_", years[y], ".png")), width = 8, height =4)
    
  }
  
  
  setwd(file.path(DIR$survey_figs, "RasterCatch"))
  # Converting .png files in one .gif image using ImageMagick
  system("magick -delay 80 *.png catch.gif")
  setwd(DIR$R)
}

## what countries were doing the surveys when
table(catch_sf$Country, catch_sf$Year)

## see if we can recreate the index data sets from raw data
catch_sf$CPUE = catch_sf$catch
table(catch_sf$Geographic.Area.Name %in% index_by_area_df$Geographic.Area.Name)
table(is.na(catch_sf$Geographic.Area.Name %in% index_by_area_df$Geographic.Area.Name))
table(index_by_area_df$Geographic.Area.Name %in% catch_sf$Geographic.Area.Name)

## attach area to catch df
catch_sf$Geographic.Area.Size..km. = index_by_area_df$Geographic.Area.Size..km.[match(catch_sf$Geographic.Area.Name, index_by_area_df$Geographic.Area.Name)]
check_nas = subset(catch_sf, subset = is.na(catch_sf$Geographic.Area.Size..km.))
table(check_nas$Geographic.Area.Name)
table(check_nas$Depth.Stratum)
##
table(index_by_area_df$Geographic.Area.Name, index_by_area_df$Geographic.Area.Size..km.)
## ask about 
## crude index
#catch_index_df = catch_sf %>% group_by(Year) %>% summarise(crude_ndx = gm_mean(CPUE))#

#png(filename = file.path(DIR$survey_figs, "quick_and_dirty_comparison.png"), width = 7, height = 5, res = 250, units = "in")
#plot(catch_index_df$Year, normalise(catch_index_df$crude_ndx), type = "o", xlab = "Year", ylab = "Normalised index (RPN)", lwd = 3)
#lines(index_df$Year[index_df$Country != "Japan"], normalise(index_df$CPUE[index_df$Country != "Japan"]), lwd = 3, col = "blue", lty = 2)
#legend('topleft', legend = c("Geometric mean", "Input index"), col= c("black","blue"), lwd = 3)
#dev.off()

## how many zeros?
table(catch_sf$catch == 0)
table(catch_sf$CPUE == 0)

################
## Characterise Age samples
###############
table(age_df_full$Country)
table(age_df_full$Year.x)

head(age_df_full)
## make it a sf object
survey_age_sf = st_as_sf(x = age_df_full,                         
                         coords = c("Start.Longitude..DD..x", "Start.Latitude..DD..x"),
                         crs = projcrs)
## transform to more convenient projection "North American Datum 1983" 
## this removes the -180 problem
survey_age_nad83_df = st_transform(survey_age_sf, crs = st_crs(longline_areas))

survey_age_nad83_df = survey_age_nad83_df %>% st_join(longline_areas)
survey_age_nad83_df$area_lab = obs_catch_labs$actual_lab[match(survey_age_nad83_df$NAME, obs_catch_labs$obs_catch_lab)]

survey_age_df_for_AF = survey_age_nad83_df  %>% st_drop_geometry()
saveRDS(survey_age_df_for_AF, file = file.path("Data", "Survey", "age_df_for_AF.RDS"))

## Plot samples of aged fish by sex and year
age_samples_by_year_and_sex = survey_age_df %>% group_by(Year, Sex.Description) %>% summarise(Age_samples = sum(!is.na(Age))) %>%
  pivot_wider(values_from = Age_samples, id_cols = Year, names_from = Sex.Description)

ggplot(age_samples_by_year_and_sex) +
  geom_line(aes(x = Year, y = Male, col = "Male"), linewidth = 1.1) +
  geom_point(aes(x = Year, y = Male, col = "Male", shape = "Male"), size = 2.1) +
  geom_line(aes(x = Year, y = female, col = "Female"), linewidth = 1.1, linetype = "dashed") +
  geom_point(aes(x = Year, y = female, col = "Female", shape = "Female"), size = 2.1) +
  labs(x = "Year", y = "Age samples", col = "Sex", shape = "Sex")
ggsave(filename = file.path(DIR$survey_figs, "age_samples_by_year_and_sex.png"), width = 8, height =4)

## Plot spatial distribution of aged fish
all_age_locations = ggplot() +
  geom_sf(data = longline_areas, aes(col = "Longline"), fill = NA, linetype = "dashed", linewidth = 1) +
  geom_sf(data = world, fill = "darkgreen") +
  geom_sf(data = survey_age_nad83_df, fill = "black") +
  coord_sf(crs  = st_crs(longline_areas), xlim = LL_xlim, ylim = LL_ylim) +
  labs(x = "", y = "", col = "") +
  theme(axis.text = element_blank(),
        legend.position = "bottom") +
  guides(col = "none")
ggsave(plot = all_age_locations, filename = file.path(DIR$survey_figs, "spatial_distribution_of_age_samples.png"), width = 8, height =4)


## Plot samples of aged fish by sex and year
age_samples_by_year_region_sex = survey_age_df %>% group_by(Year, Sex.Description, `NPFMC.Sablefish.Mgmt.Area`) %>% summarise(Age_samples = sum(!is.na(Age))) %>%
  pivot_wider(values_from = Age_samples, id_cols = c(Year,NPFMC.Sablefish.Mgmt.Area), names_from = Sex.Description)
ggplot(age_samples_by_year_region_sex) +
  geom_line(aes(x = Year, y = Male, col = "Male"), linewidth = 1.1) +
  geom_point(aes(x = Year, y = Male, col = "Male", shape = "Male"), size = 2.1) +
  geom_line(aes(x = Year, y = female, col = "Female"), linewidth = 1.1, linetype = "dashed") +
  geom_point(aes(x = Year, y = female, col = "Female", shape = "Female"), size = 2.1) +
  labs(x = "Year", y = "Age samples", col = "Sex", shape = "Sex") +
  facet_wrap(~NPFMC.Sablefish.Mgmt.Area)
ggsave(filename = file.path(DIR$survey_figs, "age_samples_by_year_region_sex.png"), width = 8, height =6)

## plot ECDF for raw Age-frequencies
age_props_by_year_and_sex = survey_age_df %>% filter(Sex != 3) %>% group_by(Year, Sex.Description) %>% mutate(age_composition = Age / sum(Age, na.rm = T))
age_props_by_year_and_sex$age_composition = age_props_by_year_and_sex$age_composition %>% replace_na(0)
age_props_by_year_and_sex = age_props_by_year_and_sex %>% group_by(Year, Sex.Description) %>% mutate(cdf = cumsum(age_composition))

## couldnt get this plot to do what I wanted, something not right with the cumsum
#ggplot(age_props_by_year_and_sex, aes(x=Age, color=factor(Sex))) +
#  geom_line(aes(y = cdf)) +
#  facet_wrap(~Year)

## Simple empircal
ggplot(survey_age_df %>% filter(Sex != 3), aes(x = Age, col = Year, group = Year)) + 
  stat_ecdf() +
  labs(x = "Age", y = "ECDF") +
  xlim(0,60) +
  facet_wrap(~Sex.Description)
ggsave(filename = file.path(DIR$survey_figs, "age_samples_by_year_CDF.png"), width = 6, height =4)

ggplot(survey_age_df %>% filter(Sex != 3), aes(x = Age, col = Sex.Description, group = Sex.Description, linetype = Sex.Description)) + 
  stat_ecdf() +
  xlim(0,60) + 
  labs(col = "Sex", linetype ="Sex", y = "Cumulative AF") +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "age_samples_by_year_sex_CDF.png"), width = 10, height =8)

## by region
ggplot(survey_age_df %>% filter(Sex != 3), aes(x = Age, col = Sex.Description, group = Sex.Description, linetype = Sex.Description)) + 
  stat_ecdf() +
  xlim(0,60) + 
  labs(col = "Sex", linetype ="Sex", y = "Cumulative AF") +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "age_samples_by_year_sex_CDF.png"), width = 10, height =8)

##
## by region
age_props_by_year_region_sex = survey_age_df %>% filter(Sex != 3) %>% group_by(Year, Sex.Description, NPFMC.Sablefish.Mgmt.Area) %>% mutate(age_composition = Age / sum(Age, na.rm = T))
age_props_by_year_region_sex$age_composition = age_props_by_year_region_sex$age_composition %>% replace_na(0)

ggplot(age_props_by_year_region_sex %>% filter(Sex.Description == "Male"), aes(x = Age, fill = NPFMC.Sablefish.Mgmt.Area, col = NPFMC.Sablefish.Mgmt.Area, y = age_composition, alpha = 0.3)) +
  geom_bar(stat = "identity") +
  labs(col = "", fill = "", y = "Raw AF") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  guides(alpha = "none") +
  xlim(0,40) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_AFs_by_year_region_male.png"), width = 10, height =10)


ggplot(age_props_by_year_region_sex %>% filter(Sex.Description == "female"), aes(x = Age, fill = NPFMC.Sablefish.Mgmt.Area, col = NPFMC.Sablefish.Mgmt.Area, y = age_composition, alpha = 0.3)) +
  geom_bar(stat = "identity") +
  labs(col = "", fill = "", y = "Raw AF") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  guides(alpha = "none") +
  xlim(0,40) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_AFs_by_year_region_female.png"), width = 10, height =10)

## repeat teh above plots but with ECDF
ggplot(survey_age_df %>% filter(Sex.Description == "female"), aes(x = Age, col = NPFMC.Sablefish.Mgmt.Area, group = NPFMC.Sablefish.Mgmt.Area)) + 
  stat_ecdf() +
  labs(x = "Age", y = "ECDF", col = "") +
  xlim(0,60) +
  theme(legend.position = "bottom")+
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_ecdf_AFs_by_year_region_female.png"), width = 10, height =10)

ggplot(survey_age_df %>% filter(Sex.Description == "Male"), aes(x = Age, col = NPFMC.Sablefish.Mgmt.Area, group = NPFMC.Sablefish.Mgmt.Area)) + 
  stat_ecdf() +
  labs(x = "Age", y = "ECDF", col = "") +
  theme(legend.position = "bottom")+
  xlim(0,60) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_ecdf_AFs_by_year_region_male.png"), width = 10, height =10)



################
## Characterise Length samples
###############
head(survey_length_df)
table(survey_length_df$Country)
## make it a sf object
survey_length_sf = st_as_sf(x = length_df_full,                         
                            coords = c("Start.Longitude..DD..x", "Start.Latitude..DD..x"),
                            crs = projcrs)
## transform to more convenient projection "North American Datum 1983" 
## this removes the -180 problem
survey_length_nad83_df = st_transform(survey_length_sf, crs = st_crs(longline_areas))

survey_length_nad83_df = survey_length_nad83_df %>% st_join(longline_areas)
survey_length_nad83_df$area_lab = obs_catch_labs$actual_lab[match(survey_length_nad83_df$NAME, obs_catch_labs$obs_catch_lab)]

## Plot samples of aged fish by sex and year
length_2cm_breaks = seq(from = min(survey_length_df$Length) - 1, to = max(survey_length_df$Length) + 1, by = 2)
survey_length_df$length_2_cm_breaks = cut(survey_length_df$Length, breaks = length_2cm_breaks)


length_samples_by_year_and_sex = survey_length_df %>% group_by(Year, Sex) %>% summarise(Length_samples = sum(Frequencey, na.rm = T)) %>%
  pivot_wider(values_from = Length_samples, id_cols = Year, names_from = Sex)

ggplot(length_samples_by_year_and_sex) +
  geom_line(aes(x = Year, y = `1`, col = "Male"), linewidth = 1.1) +
  geom_point(aes(x = Year, y = `1`, col = "Male", shape = "Male"), size = 2.1) +
  geom_line(aes(x = Year, y = `2`, col = "Female"), linewidth = 1.1, linetype = "dashed") +
  geom_point(aes(x = Year, y = `2`, col = "Female", shape = "Female"), size = 2.1) +
  geom_line(aes(x = Year, y = `3`, col = "Unknown"), linewidth = 1.1, linetype = "dotted") +
  geom_point(aes(x = Year, y = `3`, col = "Unknown", shape = "Unknown"), size = 2.1) +
  labs(x = "Year", y = "Length samples", col = "Sex", shape = "Sex")
ggsave(filename = file.path(DIR$survey_figs, "length_samples_by_year_and_sex.png"), width = 8, height =4)


## Plot samples of lengthed fish by sex and year
length_samples_by_year_region_sex = survey_length_df %>% group_by(Year, Sex, `NPFMC.Sablefish.Management.Area`) %>% summarise(Length_samples = sum(Frequencey, na.rm = T)) %>%
  pivot_wider(values_from = Length_samples, id_cols = c(Year,NPFMC.Sablefish.Management.Area), names_from = Sex)
ggplot(length_samples_by_year_region_sex) +
  geom_line(aes(x = Year, y = `1`, col = "Male"), linewidth = 1.1) +
  geom_point(aes(x = Year, y = `1`, col = "Male", shape = "Male"), size = 2.1) +
  geom_line(aes(x = Year, y = `2`, col = "Female"), linewidth = 1.1, linetype = "dashed") +
  geom_point(aes(x = Year, y = `2`, col = "Female", shape = "Female"), size = 2.1) +
  labs(x = "Year", y = "Length samples", col = "Sex", shape = "Sex") +
  facet_wrap(~NPFMC.Sablefish.Management.Area)
ggsave(filename = file.path(DIR$survey_figs, "length_samples_by_year_region_sex.png"), width = 8, height =6)


## plot ECDF for raw Length-frequencies
length_pooled_by_year_and_sex = survey_length_df %>% filter(Sex != 3) %>% group_by(Year, Sex, Length) %>% summarise(pooled_frequency = sum(Frequencey, na.rm = T))
length_props_by_year_and_sex = length_pooled_by_year_and_sex %>% group_by(Year, Sex) %>% mutate(length_composition = pooled_frequency / sum(pooled_frequency))
length_props_by_year_and_sex$length_composition = length_props_by_year_and_sex$length_composition %>% replace_na(0)
length_props_by_year_and_sex = length_props_by_year_and_sex %>% group_by(Year, Sex) %>% mutate(cdf = cumsum(length_composition))


## by region
length_pooled_by_year_region_sex = survey_length_df %>% filter(Sex != 3) %>% group_by(Year, Sex, NPFMC.Sablefish.Management.Area, Length) %>% summarise(pooled_frequency = sum(Frequencey, na.rm = T))
length_props_by_year_region_sex = length_pooled_by_year_region_sex %>%  group_by(Year, Sex, NPFMC.Sablefish.Management.Area) %>% mutate(length_composition = pooled_frequency / sum(pooled_frequency), cdf_LF = cumsum(pooled_frequency) / sum(pooled_frequency))
length_props_by_year_region_sex$length_composition = length_props_by_year_region_sex$length_composition %>% replace_na(0)

ggplot(length_props_by_year_region_sex %>% filter(Sex == 1), aes(x = Length, fill = NPFMC.Sablefish.Management.Area, col = NPFMC.Sablefish.Management.Area, y = length_composition, alpha = 0.3)) +
  geom_bar(stat = "identity") +
  labs(col = "", fill = "", y = "Raw LF") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  guides(alpha = "none") +
  xlim(25,100) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_LFs_by_year_region_male.png"), width = 10, height =10)


ggplot(length_props_by_year_region_sex %>% filter(Sex == 2), aes(x = Length, fill = NPFMC.Sablefish.Management.Area, col = NPFMC.Sablefish.Management.Area, y = length_composition, alpha = 0.3)) +
  geom_bar(stat = "identity") +
  labs(col = "", fill = "", y = "Raw LF") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  guides(alpha = "none") +
  xlim(25,100) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_LFs_by_year_region_female.png"), width = 10, height =10)

## repeat teh above plots but with ECDF
ggplot(length_props_by_year_region_sex %>% filter(Sex == 2), aes(x = Length, col = NPFMC.Sablefish.Management.Area, group = NPFMC.Sablefish.Management.Area)) + 
  geom_line(aes(y = cdf_LF)) +
  labs(x = "Length", y = "ECDF", col = "") +
  xlim(40,90) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_ecdf_LFs_by_year_region_female.png"), width = 10, height =10)

ggplot(length_props_by_year_region_sex %>% filter(Sex == 1), aes(x = Length, col = NPFMC.Sablefish.Management.Area, group = NPFMC.Sablefish.Management.Area)) + 
  geom_line(aes(y = cdf_LF)) +
  labs(x = "Length", y = "ECDF", col = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14)
  ) +
  xlim(40,90) +
  facet_wrap(~Year)
ggsave(filename = file.path(DIR$survey_figs, "raw_ecdf_LFs_by_year_region_male.png"), width = 10, height =10)

##########
## Calculate scaled LFs
##########
length_df_for_LF = survey_length_nad83_df  %>% st_drop_geometry()
length_df_for_LF = length_df_for_LF %>% group_by(event_id) %>% mutate(sampling_prop = sum(Frequencey) / unique(catch))
sum_df = length_df_for_LF %>% group_by(event_id) %>% summarise(LF_samples = sum(Frequencey), event_catch =  mean(catch), sampling_prop = LF_samples/event_catch)
## no catch recorded for these stations?
na_sum_df = sum_df %>% filter(is.na(sampling_prop))
station_catch %>% filter(event_id %in% na_sum_df$event_id)
##
catch_df  %>% filter(event_id %in% na_sum_df$event_id)
## Ask someone about this. Onwards and upwards
hist(length_df_for_LF$sampling_prop)
summary(length_df_for_LF$sampling_prop)

## look at events where sampled more LFs than in catch
print(sum_df %>% filter(sampling_prop > 1), n = 140)
## most of these are old < 1998 so I am going to set the sampling Prop = 1 
## This assumes catch = LF samples
length_df_for_LF = length_df_for_LF %>% mutate(ifelse(sampling_prop > 1, 1, sampling_prop))

## Scale LFs to event level
length_df_for_LF$event_lvl_frequency = length_df_for_LF$Frequencey / length_df_for_LF$sampling_prop
hist(length_df_for_LF$sampling_prop, breaks = 40, xlim = c(0,6))
saveRDS(length_df_for_LF, file = file.path("Data", "Survey", "length_df_for_LF.RDS"))

## take the mean LF within the region of interest
LF_by_area_Year_sex = length_df_for_LF %>% group_by(area_lab, Year.x, Sex, Length) %>% summarise(mean_length = sum(event_lvl_frequency) / length(unique(event_id)))
LF_by_area_Year_sex = LF_by_area_Year_sex %>% group_by(area_lab, Year.x, Sex) %>% mutate(prop = mean_length / sum(mean_length)) 
head(survey_length_df)
colnames(survey_length_df)
print(length_df_for_LF %>% group_by(area_lab, Year.x) %>% summarise(sample_size = sum(Frequencey)) %>% pivot_wider(id_cols = Year.x, names_from = area_lab, values_from = sample_size), n = 44)
## visualize
ggplot(LF_by_area_Year_sex %>% group_by(area_lab, Year.x, Sex) %>%  
         arrange((Length)) %>%
         mutate(cumsum=cumsum(prop)) %>% filter(Year.x %in% 1979:2008, Sex == 1), aes(x = Length, y = (cumsum), col = area_lab, linetype = area_lab)) +
  geom_line(linewidth = 1.1) +  
  labs(x = "Length (cm)") +
  facet_wrap(~Year.x, ncol = 4) +
  ggtitle("Male LL LF") +
  xlim(40, 80) +
  theme_bw()
ggplot(LF_by_area_Year_sex %>% group_by(area_lab, Year.x, Sex) %>%  
         arrange((Length)) %>%
         mutate(cumsum=cumsum(prop)) %>% filter(Year.x %in% 1979:2008, Sex == 2), aes(x = Length, y = (cumsum), col = area_lab, linetype = area_lab)) +
  geom_line(linewidth = 1.1) +  
  labs(x = "Length (cm)") +
  facet_wrap(~Year.x, ncol = 4) +
  ggtitle("Female LL LF") +
  xlim(40, 80) +
  theme_bw()

## over both sexes
LF_by_area_Year_sex = length_df_for_LF %>% group_by(area_lab, Year.x, Length) %>% summarise(mean_length = sum(event_lvl_frequency) / length(unique(event_id)))
LF_by_area_Year_sex = LF_by_area_Year_sex %>% group_by(area_lab, Year.x) %>% mutate(prop = mean_length / sum(mean_length)) 
head(survey_length_df)
colnames(survey_length_df)
print(length_df_for_LF %>% group_by(area_lab, Year.x) %>% summarise(sample_size = sum(Frequencey)) %>% pivot_wider(id_cols = Year.x, names_from = area_lab, values_from = sample_size), n = 44)
## visualize
ggplot(LF_by_area_Year_sex %>% group_by(area_lab, Year.x) %>%  
         arrange((Length)) %>%
         mutate(cumsum=cumsum(prop)) %>% filter(Year.x %in% 1979:2002), aes(x = Length, y = (cumsum), col = area_lab, linetype = area_lab)) +
  geom_line(linewidth = 1.1) +  
  labs(x = "Length (cm)", y = "Cumulative proportion", col = "Area", linetype = "Area") +
  facet_wrap(~Year.x, ncol = 4) +
  ggtitle("Sex aggregated LF") +
  xlim(40, 80) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 12))
ggsave(filename = file.path(DIR$survey_figs, "scale_LF_sex_agg.png"), width = 10, height =12)
## visualize
ggplot(LF_by_area_Year_sex %>% group_by(area_lab, Year.x) %>%  
         arrange((Length)) %>%
         mutate(cumsum=cumsum(prop)) %>% filter(Year.x %in% 2003:2022), aes(x = Length, y = (cumsum), col = area_lab, linetype = area_lab)) +
  geom_line(linewidth = 1.1) +  
  labs(x = "Length (cm)", y = "Cumulative proportion", col = "Area", linetype = "Area") +
  facet_wrap(~Year.x, ncol = 4) +
  ggtitle("Sex aggregated LF") +
  xlim(40, 80) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 12))
ggsave(filename = file.path(DIR$survey_figs, "scale_LF_sex_agg_2.png"), width = 10, height =12)

################
## Growth
###############
## look for errors
## number of aged fish by sex
age_samples = survey_age_df %>% group_by(Year) %>% summarise(n_age_samples = sum(!is.na(Age)))
print(age_samples, n = 32)
decades = seq(from = 1980, to = 2030, by = 10)
decade_label = c(paste0(decades[1:(length(decades) - 1)], "-",paste0(decades[2:length(decades)])))
survey_age_df$decade = cut(survey_age_df$Year, breaks = decades, labels = decade_label)


ggplot(survey_age_df, aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = Sex.Description, shape = Sex.Description)) +
  guides(alpha = "none")
## look at suspicous 0 aged fish
ggplot(survey_age_df %>% filter(Age == 0), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = Sex.Description, shape = Sex.Description)) +
  guides(alpha = "none")

## re do dropping 0 year old fish
plt = ggplot(survey_age_df %>% filter(Age > 0), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = Sex.Description, shape = Sex.Description)) +
  labs(x = "Age", y = "Length (cm)", col = "Sex", shape = "Sex", title = "") +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none") +
  theme(legend.position = "bottom")
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_growth_by_sex.png"), width = 6, height =6)

## Male by area
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "Male"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = NPFMC.Sablefish.Mgmt.Area, shape = NPFMC.Sablefish.Mgmt.Area)) +
  labs(x = "Age", y = "Length (cm)", col = "Region", shape = "Region", title = "Male") +
  geom_smooth(aes(group = NPFMC.Sablefish.Mgmt.Area, col = NPFMC.Sablefish.Mgmt.Area, linetype = NPFMC.Sablefish.Mgmt.Area), se = FALSE) +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none",linetype = 'none') +
  theme(legend.position = "bottom") 
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_male_growth_by_region.png"), width = 6, height =6)

## Female by area
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "female"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = NPFMC.Sablefish.Mgmt.Area, shape = NPFMC.Sablefish.Mgmt.Area)) +
  labs(x = "Age", y = "Length (cm)", col = "Region", shape = "Region", title = "Female") +
  xlim(0, 75) +
  geom_smooth(aes(group = NPFMC.Sablefish.Mgmt.Area, col = NPFMC.Sablefish.Mgmt.Area, linetype = NPFMC.Sablefish.Mgmt.Area), se = FALSE) +
  ylim(35, 120) +
  guides(alpha = "none",linetype = 'none') +
  theme(legend.position = "bottom") 
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_female_growth_by_region.png"), width = 6, height =6)

## Male by Year
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "Male"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = Year)) +
  labs(x = "Age", y = "Length (cm)", col = "Year") +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none") +
  theme(legend.position = "bottom")
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_male_growth_by_year.png"), width = 6, height =6)

## female by Year
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "female"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = Year)) +
  labs(x = "Age", y = "Length (cm)", col = "Year") +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none") +
  theme(legend.position = "bottom")
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_female_growth_by_year.png"), width = 6, height =6)

## Male by decade
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "Male"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = decade)) +
  labs(x = "Age", y = "Length (cm)", col = "Decade") +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none") +
  theme(legend.position = "bottom")
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_male_growth_by_decade.png"), width = 6, height =6)

## female by decade
plt = ggplot(survey_age_df %>% filter(Age > 0, Sex.Description == "female"), aes(x = Age, y = Length..cm., alpha = 0.5)) +
  geom_point(aes(col = decade)) +
  labs(x = "Age", y = "Length (cm)", col = "Decade") +
  xlim(0, 75) +
  ylim(35, 120) +
  guides(alpha = "none") +
  theme(legend.position = "bottom")
ggsave(plot = plt, filename = file.path(DIR$survey_figs, "raw_female_growth_by_decade.png"), width = 6, height =6)


###############
## Explore spatial structure in 
## Survey length data.
###############
head(length_df_for_LF)

## look at the difference in time
length_df_for_LF$Haul.Date.x = as.Date(length_df_for_LF$Haul.Date.x, format = "%d-%h-%y")
length_df_for_LF$LonDD.End_wrap180 = west2east(length_df_for_LF$End.Longitude..DD..x)

## reformat length data to be compatible with Haikun's library
length_df_for_LF$quarter = quarters(length_df_for_LF$Haul.Date.x)
length_df_for_LF$month = month(length_df_for_LF$Haul.Date.x)

## aggregate LF data at some spatial resolution
x_res = 2 # degrees
y_res = 2 # degrees

lon_ext = c(170, 240)
lat_ext = c(40, 60)
#range(observer_length_df_full$LonDD.End_wrap180, na.rm = T)
#range(observer_length_df_full$LatDD.End, na.rm = T)

LF_grd = raster(vals = NA, resolution = c(x_res, y_res), ext = extent(c(lon_ext, lat_ext)))


years = sort(unique(length_df_for_LF$Year.x))
seasons = unique(length_df_for_LF$quarter)
len_range = range(length_df_for_LF$Length, na.rm = T)
LF_res = 2 ## may need to explore with this.
LF_min_class = seq(from = 40, to = 80, by = LF_res) ## make the last length bin aplus group
LF_max_class = seq(from = 39.99, to = 80, by = LF_res) ## make the last length bin aplus group
LF_min_class = c(0, LF_min_class)
LF_max_class = c(LF_max_class, 200)
length(LF_min_class)== length(LF_max_class)
length_labels = paste0(LF_min_class, "-", c(LF_min_class[2:length(LF_min_class)], "80+"))
## reformat data
sample_size = matrix(NA, nrow = length(years), ncol = length(seasons), dimnames = list(years, seasons))
#gear_type = "LONGLINER"
if(FALSE) {
  full_df = NULL
  for(y_ndx in 1:length(years)) {
    for(s_ndx in 1:length(seasons)) {
      this_df = length_df_for_LF %>% filter(Year.x == years[y_ndx], quarter == seasons[s_ndx], !is.na(LonDD.End_wrap180), !is.na(End.Latitude..DD..x))
      if(!nrow(this_df) > 0)
        next;
      sample_size[y_ndx, s_ndx] = sum(this_df$Frequencey, na.rm = T)
      ## rasterise for each length bin! quick and dirty
      for(len_ndx in 1:length(LF_min_class)) {
        this_length_df = this_df %>% filter(Length > LF_min_class[len_ndx], Length < LF_max_class[len_ndx])
        if(nrow(this_length_df) > 0) {
          this_rster  = rasterize(x = this_length_df[, c("LonDD.End_wrap180","End.Latitude..DD..x")], y = LF_grd, field = this_length_df$Frequencey, fun = sum, na.rm = T)
          this_rster_df = as.data.frame(this_rster, xy = TRUE)
          this_rster_df = subset(this_rster_df, subset = !is.na(this_rster_df$layer))
          this_rster_df$length_labels = length_labels[len_ndx]
          this_rster_df$year = years[y_ndx]
          this_rster_df$season = seasons[s_ndx]
          full_df = rbind(full_df, this_rster_df)
        }
      }
    }
  }
  ## save this
  #saveRDS(object = full_df, file = file.path("Data", "survey_LF_reformatted_2_res.RDS"))
}
## read in
full_df = readRDS(file = file.path("Data", "survey_LF_reformatted_2_res.RDS"))
full_df$quarter = as.numeric(substring(full_df$season, first = 2))
full_df$lat = full_df$y
full_df$lon = full_df$x

## covert numbers to proportions
full_df= full_df %>% group_by(year, quarter, lat, lon) %>% mutate(sample_size =sum(layer, na.rm = T), length_props = layer / sample_size)
# check they sum = 1
test = full_df %>% group_by(year, quarter, lat, lon) %>% summarise(check = sum(length_props))
table(test$check == 1)
## pivot wider so its in correct format
LF_long = full_df %>% pivot_wider(id_cols = c(year, quarter, lat, lon, sample_size), names_from = length_labels, values_from = length_props, values_fill = 0)

##
fcol = 6
lcol = 26
bins = LF_min_class[-1]
Nsplit = 3
save_dir = file.path(DIR$survey_figs, "LF_length_analysis")

## check rows sum to one
row_sum <- apply(LF_long[,fcol:lcol],1,sum)
bad_ndx = which(abs(row_sum-1)>0.05)
length(bad_ndx)
LF_long[bad_ndx, ]
row_sum[bad_ndx]
## drop them
LF_long = subset(LF_long, subset = abs(row_sum-1)<=0.05)

make.lf.map(LF = LF_long, fcol, lcol, bins, save_dir, plot_name = "make_lf_map.png", plot_format = "png")
ggsave(filename = file.path(DIR$survey_figs, "LF_length_analysis", "make_lf_map.png"), height = 7, width = 12)
make.meanl.map(LF_long, fcol, lcol, bins, save_dir, s =13)
ggsave(filename = file.path(DIR$survey_figs, "LF_length_analysis", "make.meanl.map.png"), height = 7, width = 12)

## run regression tree
## turn off quarter and year
LF_tree = run_regression_tree(LF = LF_long, fcol, lcol, bins,Nsplit =2, save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=F,year=F,include_dummy=FALSE,pdf=FALSE)
#
LF_tree$Record
make.split.map(LF = LF_tree$LF,Nsplit =2,save_dir = save_dir, s = 2)
## turn on year
LF_tree_yr = run_regression_tree(LF = LF_long, fcol, lcol, bins,Nsplit =5, save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=F,year=T,include_dummy=FALSE,pdf=FALSE)
LF_tree_yr$Record
make.split.map(LF = LF_tree_yr$LF,Nsplit =4,save_dir = save_dir, s = 2)

## turn off year and quarter on
LF_tree_qtr = run_regression_tree(LF = LF_long, fcol, lcol, bins,Nsplit =5, save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=T,year=F,include_dummy=FALSE,pdf=FALSE)
LF_tree_qtr$Record
make.split.map(LF = LF_tree_qtr$LF,Nsplit =4,save_dir = save_dir, s = 2)

## turn on year and quarter
LF_tree_yr_qtr = run_regression_tree(LF = LF_long, fcol, lcol, bins,Nsplit =5, save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=T,year=T,include_dummy=FALSE,pdf=FALSE)
LF_tree_yr_qtr$Record
make.split.map(LF = LF_tree_yr_qtr$LF,Nsplit =4,save_dir = save_dir, s = 2)

saveRDS(LF_tree_yr_qtr$Record, file = file.path(DIR$book, "Data", "survey_LF_split_algorithm.RDS"))

split_coords_1 = matrix(c(206, 60, 206, 40), byrow = T, ncol = 2)
split_coords_2 = matrix(c(178, 60, 178, 40), byrow = T, ncol = 2)
split_coords_3 = matrix(c(206, 60, 206, 40), byrow = T, ncol = 2)
split_lons_1 = st_cast(st_sfc(st_linestring(split_coords_1)), "LINESTRING")
split_lons_1 = st_set_crs(split_lons_1, (projcrs))
split_lons_3 = st_cast(st_sfc(st_linestring(split_coords_3)), "LINESTRING")
split_lons_3 = st_set_crs(split_lons_3, (projcrs))
split_lons_2 = st_cast(st_sfc(st_linestring(split_coords_2)), "LINESTRING")
split_lons_2 = st_set_crs(split_lons_2, (projcrs))

## visualise these splits on a map.
all_station_plot = ggplot() +
  # geom_sf(data = longline_areas, aes(col = "Boundarys"), fill = NA, linetype = "dashed", linewidth = 1) +
  # geom_sf(data = world, fill = "darkgreen") +
  geom_sf(data = split_lons_1, aes(col = "Split: 2"), linetype = "solid", linewidth = 1) +
  geom_sf(data = split_lons_2, aes(col = "Split: 3"), linetype = "solid", linewidth = 1) +
  geom_sf(data = split_lons_3, aes(col = "Split: 5"), linetype = "dashed", linewidth = 1) +
  coord_sf(crs  = st_crs(longline_areas), xlim = LL_xlim, ylim = LL_ylim)  +
  labs(x = "", y = "", col = "") +
  theme(axis.text = element_blank(),
        legend.position = "bottom") 

all_station_plot
ggsave(plot = all_station_plot, filename = file.path(DIR$survey_figs, "survey_split_analysis.png"), width = 8, height =4)
