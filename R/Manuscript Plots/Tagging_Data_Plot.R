# Purpose: To visualize tagging data
# Creator: Matthew LH. Cheng (UAF - CFOS)
# date 9/23/24


# Set up ------------------------------------------------------------------
library(tidyverse)
library(here)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(ggpubr)

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
    round(Time_Liberty) == 0 & round(Time_Liberty) < 6 ~ "0 - 5 Years",
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

# Plot --------------------------------------------------------------------

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
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.65, 0.85), legend.background = element_blank(), legend.box = element_blank(),
        legend.spacing.x = unit(1, "cm")) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Recovery Year", y = "Number of Recoveries", color = "Gear Type", lty = "Gear Type")

# Time at liberty
liberty_hist = ggplot(recovery_plot_df, aes(x = Time_Liberty)) +
  geom_histogram(bins = 40, color = 'black', alpha = 0.5) +
  theme_bw(base_size = 10) +
  labs(x = "Time at Liberty (Years)", y = "Number of Recoveries")
  
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
  labs(x = "Longitude", y = "Latitude", fill = "Alaska Sablefish Management Regions", color = "Tag Recovery Direction")

# Output plot
legend = get_legend(liberty_map_plot + theme(legend.position="right")) # get legend from map
top = plot_grid(tag_recovery_ts, liberty_hist, ncol = 1, labels = c("A", "B")) # get plot on top
combined = plot_grid(top, liberty_map_plot, ncol = 2, rel_widths = c(0.35, 0.65), labels = c("", "C"), label_x = 0.095)
ggsave(here("Figs", "Manuscript_Plots", "Tag_Recovery_Summary.png"),
       plot_grid(combined, legend_b, ncol = 2, rel_widths = c(0.75, 0.25)))

