# To plot map and tagging data inputs 
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date 9/23/24

# Set Up ------------------------------------------------------------------
library(tidyverse)
library(here)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)

# Get data
haul_df = read.csv(file = file.path(here("data", "Tagging"), "HAUL_VIEW.csv"))
release_df = read.csv(file = file.path(here("data", "Tagging"), "RELEASE_MVIEW_Sablefish.csv"))
recovery_df = read.csv(file = file.path(here("data", "Tagging"), "RECOVERY_MVIEW_Sablefish_Groomed for Craig.csv"))

# Read in stat areas
nmfs_areas = read_sf(dsn = here("data", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas = st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas = nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas = st_shift_longitude(nmfs_areas) # shift longitude for plotting
nmfs_areas = nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A"))

# Get map for plotting
west = ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west = st_shift_longitude(west) # shift ongitude for plotting


# Data Munging ------------------------------------------------------------

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

# Define east and west general direction
recovery_df = recovery_df %>% 
  filter(TIME_OUT > 0) %>% 
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

# Plot --------------------------------------------------------------------

map_plot = ggplot() +
geom_sf(data = nmfs_areas, aes(fill = NAME, lty = GEN_NAME), alpha = 0.5, lwd = 0.2) + # Sablefish Stat Area
geom_segment(recovery_df, mapping = aes(x = HLNG2, HLAT, xend = RLNG2, yend = RLAT), alpha = 0.055, lineend = "round") + # Direction
geom_point(recovery_df, mapping = aes(x = HLNG2, y = HLAT), pch = 20, position = position_jitter(width = 0.05), alpha = 0.3, size = 2) + # Release
geom_point(recovery_df, mapping = aes(x = RLNG2, y = RLAT, color = Direction), size = 2, alpha = 0.3, pch = 18, position = position_jitter(width = 0.05)) + # Recovery 
geom_point(mean_sd_centroid, mapping = aes(x = mean_HLNG, y = mean_HLAT), pch = 19, size = 5, color = 'white') + # Release centroid
geom_errorbar(mean_sd_centroid, mapping = aes(x = mean_HLNG, y = mean_HLAT,  # Release SD
                                              ymin = mean_HLAT - SD_HLAT, ymax = mean_HLAT + SD_HLAT), lty = 2, color = 'white', lwd = 0.55) +
geom_errorbarh(mean_sd_centroid, mapping = aes(y = mean_HLAT, xmin = mean_HLNG - SD_HLNG, # Release SD
                                               xmax = mean_HLNG + SD_HLNG), lty = 2, color = 'white', lwd = 0.55) +
# geom_point(mean_sd_centroid, mapping = aes(x = mean_RLNG, y = mean_RLAT), pch = 18, size = 5, color = 'white') + # Recovery centroid
# geom_errorbar(mean_sd_centroid, mapping = aes(x = mean_RLNG, y = mean_RLAT, # Recovery SD
#                                               ymin = mean_RLAT - SD_RLAT, ymax = mean_RLAT + SD_RLAT), lty = 2, color = 'white', lwd = 0.55) +
# geom_errorbarh(mean_sd_centroid, mapping = aes(y = mean_RLAT, xmin = mean_RLNG - SD_RLNG, # Recovery SD
#                                                xmax = mean_RLNG + SD_RLNG), lty = 2, color = 'white', lwd = 0.55) +
geom_sf(data = west, lwd = 0.2, color = 'black') + # World Map
  
# Annotate with labels
annotate("text", label = "Alaska", x = 208, y = 65, size = 8) + # Alaska label
annotate("text", label = "Canada", x = 224.65, y = 63, size = 8) + # Canada label
annotate("text", label = "Russia", x = 172, y = 67, size = 8) + # Canada label
  
coord_sf(ylim = c(48, 70.5), xlim = c(165, 227)) + # Restrict Map Area
scale_color_manual(values = c("red", "blue")) +
scale_fill_manual(values = c("black", "#56B4E9", "#009E73","#F0E442", "#D55E00", "#E69F00")) +
scale_linetype_discrete(guide = "none") +
guides(color = guide_legend(override.aes = list(alpha = 1, size = 5, pch = c(19, 19))),
       fill = guide_legend (position = "bottom")) +
theme_bw(base_size = 18) +
theme(legend.position = c(0.85, 0.075), legend.box = "vertical", legend.background = element_blank()) +
labs(x = "Longitude", y = "Latitude", fill = "Alaska Sablefish Management Regions", color = "Tag Recovery Direction")

ggsave(here("figs", "Manuscript_Plots", "Map_Tag_Plot.png"),
       map_plot, width = 15)
