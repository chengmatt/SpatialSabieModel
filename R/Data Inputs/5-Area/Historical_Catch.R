#'
#'
#' Data from 1977-1989 was supplied by Kari Fenske who explored this in her Thesis
#'
#' The aim of this script is to come up with our best guess on historic catch from 1960-1977
#'
#'

# source("(00) Init.R")
library(dplyr)
library(ggplot2)
library(tidyr)
library(here)

## Read in Catch
catch_post_90 = read.csv(file = file.path("Data", "Raw", "raw_catch_file.csv"))
catch_90 <- read.csv(file.path("Data", "Inputs","1990 catch.csv"),header=T)
historic_catch = read.csv(file = file.path("Data", "Catch", "historic_split_catch.csv"))
fixed_gear_by_year = read.csv(file = file.path("Data", "Raw", "fsh1_catch_data.csv"))
trwl_gear_by_year = read.csv(file = file.path("Data", "Raw", "fsh2_catch_data.csv"))
fixed_gear_by_year$fishery = "Fixed"
trwl_gear_by_year$fishery = "Trawl"
total_by_year = trwl_gear_by_year
total_by_year$fishery = "Total"
total_by_year$catch = total_by_year$catch + fixed_gear_by_year$catch
full_gear_by_year = rbind(fixed_gear_by_year, trwl_gear_by_year, total_by_year)

ggplot(full_gear_by_year, aes(x = year, y = catch, col = fishery, linetype = fishery)) +
  geom_line(linewidth = 1.5) +
  geom_vline(xintercept = 1976.5, col = "gray60", linetype = "dashed", linewidth = 2) +
  theme_bw()


table(catch_post_90$fmp_gear)
table(catch_post_90$area)


## Reformat post 90's data
catch_post_90 = catch_post_90 %>% mutate(area = case_when(area == "EG"  ~ "EGOA",
                                                          area   == "CG"  ~ "CGOA",
                                                          area   == "WG"  ~ "WGOA",
                                                          TRUE ~ area),
                                         fmp_gear = case_when(fmp_gear == "POT" ~ "HAL",
                                                              TRUE ~ fmp_gear))
## reformat 1990 data
catch_90_lng = catch_90 %>% group_by(area, type) %>% summarise(catch = sum(catch))
catch_90_lng = catch_90_lng %>% mutate(area = case_when(area == "CG"  ~ "CGOA",
                                                        area == "WG"  ~ "WGOA",
                                                        area == "EG"  ~ "EGOA",
                                                        TRUE ~ area),
                                       fmp_gear = case_when(type == "trawl"  ~ "TRW",
                                                            type == "fixed"  ~ "HAL"))  %>%  select(!type)
catch_90_lng$year = 1990
catch_post_90_lng = catch_post_90 %>% group_by(year, area, fmp_gear) %>% summarise(catch = sum(weight_posted))
## reformat 1977-1989 data
## drop sum column
historic_catch = historic_catch %>% dplyr::select(!sum)
## join historic catch andto full catch
hist_catch_lng = historic_catch %>% pivot_longer(!Year)
hist_catch_lng$area = Reduce(c, lapply(strsplit(hist_catch_lng$name, split = "\\."), FUN = function(x) {x[2]}))
hist_catch_lng$gear = Reduce(c, lapply(strsplit(hist_catch_lng$name, split = "\\."), FUN = function(x) {x[1]}))
hist_catch_lng = hist_catch_lng %>% mutate(fmp_gear = case_when(gear == "trawl"  ~ "TRW",
                                                                gear == "fixed" ~ "HAL",
                                                                gear == "Fixed" ~ "HAL")) %>%
  rename(year = Year, catch = value) %>% ungroup() %>% dplyr::select(!c(name, gear))
hist_catch_lng$catch = hist_catch_lng$catch * 1000 ## its in kilotonnes, but because full_catc_df is in tonnes. we divide by 1000  later on

## join them all
full_catch_df = rbind(hist_catch_lng, catch_90_lng, catch_post_90_lng)
## Combine jig and OTH into HAL
full_catch_df = full_catch_df %>% mutate(fmp_gear = case_when(fmp_gear == "HAL" ~ "HAL",
                                                              fmp_gear == "TRW" ~ "TRW",
                                                              TRUE ~ "HAL")) %>%
  group_by(fmp_gear, year, area) %>% summarise(catch = sum(catch))



## Visualize
ggplot(full_catch_df) +
  geom_line(aes(x = year, y = catch, col = fmp_gear, linetype = fmp_gear), linewidth = 1.1) +
  facet_wrap(~area)


## historic proportions
all_data_props = full_catch_df %>% group_by(area, fmp_gear) %>% summarise(total_catch = sum(catch)) %>%
  group_by(fmp_gear) %>% mutate(proportion_catch_by_region = total_catch / sum(total_catch))

## use-data between 1977-1981
first_five_year_props = full_catch_df %>% filter(year < 1982) %>% group_by(area, fmp_gear) %>% summarise(total_catch = sum(catch)) %>%
  group_by(fmp_gear) %>% mutate(proportion_catch_by_region = total_catch / sum(total_catch))

## check they sum = 1 over the regions for each gear
all_data_props %>% group_by(fmp_gear) %>% summarise(sum(proportion_catch_by_region))
first_five_year_props %>% group_by(fmp_gear) %>% summarise(sum(proportion_catch_by_region))


fixed_gear_by_year$fmp_gear = "HAL"
trwl_gear_by_year$fmp_gear = "TRW"


fixed_gear_with_imputations_S1 = fixed_gear_by_year %>% inner_join(all_data_props)
trawl_gear_with_imputations_S1 = trwl_gear_by_year %>% inner_join(all_data_props)
fixed_gear_with_imputations_S1$imputed_catch = fixed_gear_with_imputations_S1$catch * fixed_gear_with_imputations_S1$proportion_catch_by_region
trawl_gear_with_imputations_S1$imputed_catch = trawl_gear_with_imputations_S1$catch * trawl_gear_with_imputations_S1$proportion_catch_by_region
fixed_gear_with_imputations_S2 = fixed_gear_by_year %>% inner_join(first_five_year_props)
trawl_gear_with_imputations_S2 = trwl_gear_by_year %>% inner_join(first_five_year_props)
fixed_gear_with_imputations_S2$imputed_catch = fixed_gear_with_imputations_S2$catch * fixed_gear_with_imputations_S2$proportion_catch_by_region
trawl_gear_with_imputations_S2$imputed_catch = trawl_gear_with_imputations_S2$catch * trawl_gear_with_imputations_S2$proportion_catch_by_region


ggplot() +
  geom_line(data = fixed_gear_with_imputations_S1, aes(x = year, y = imputed_catch * 1000, col = "OM1", linetype = "OM1"), linewidth = 1.2) +
  geom_line(data = fixed_gear_with_imputations_S2, aes(x = year, y = imputed_catch * 1000, col = "OM2", linetype = "OM2"), linewidth = 1.2) +
  geom_line(data = full_catch_df %>% filter(fmp_gear == "HAL"), aes(x = year, y = catch, col = "reported", linetype = "reported"), linewidth = 1.2) +
  facet_wrap(~area) +
  labs(x = "Year", y = "Catch (t)", col = "", linetype = "", title= "Fixed fishery")+
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave(filename = file.path("Output", "5-Area", "Fixed_gear_catch_area_imputation_bzero_simulation.png"), width = 10, height =7)


ggplot() +
  geom_line(data = trawl_gear_with_imputations_S1, aes(x = year, y = imputed_catch * 1000, col = "OM1", linetype = "OM1"), linewidth = 1.2) +
  geom_line(data = trawl_gear_with_imputations_S2, aes(x = year, y = imputed_catch * 1000, col = "OM2", linetype = "OM2"), linewidth = 1.2) +
  geom_line(data = full_catch_df %>% filter(fmp_gear == "TRW"), aes(x = year, y = catch, col = "reported", linetype = "reported"), linewidth = 1.2) +
  facet_wrap(~area) +
  labs(x = "Year", y = "Catch (t)", col = "", linetype = "", title= "Trawl fishery")+
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave(filename = file.path("Output", "5-Area", "Trawl_gear_catch_area_imputation_bzero_simulation.png"), width = 10, height =7)

## check they sum to the correct year values
ggplot() +
  geom_line(data = fixed_gear_with_imputations_S1 %>% group_by(year) %>% summarise(yearly_catch = sum(imputed_catch * 1000)), aes(x = year, y = yearly_catch, col = "OM1", linetype = "OM1"), linewidth = 1.2) +
  geom_line(data = fixed_gear_with_imputations_S2 %>% group_by(year) %>% summarise(yearly_catch = sum(imputed_catch * 1000)), aes(x = year, y = yearly_catch, col = "OM2", linetype = "OM2"), linewidth = 1.2) +
  geom_line(data = fixed_gear_by_year, aes(x = year, y = catch * 1000, col = "reported", linetype = "reported"), linewidth = 1.2) +
  #geom_line(data = full_catch_df %>% filter(fmp_gear == "HAL") %>% group_by(year) %>% summarise(yearly_catch = sum(catch)), aes(x = year, y = yearly_catch, col = "reported", linetype = "reported"), linewidth = 1.2) +
  labs(x = "Year", y = "Catch", col = "", linetype = "") +
  theme_bw()
ggplot() +
  geom_line(data = trawl_gear_with_imputations_S1 %>% group_by(year) %>% summarise(yearly_catch = sum(imputed_catch * 1000)), aes(x = year, y = yearly_catch, col = "OM1", linetype = "OM1"), linewidth = 1.2) +
  geom_line(data = trawl_gear_with_imputations_S2 %>% group_by(year) %>% summarise(yearly_catch = sum(imputed_catch * 1000)), aes(x = year, y = yearly_catch, col = "OM2", linetype = "OM2"), linewidth = 1.2) +
  geom_line(data = trwl_gear_by_year, aes(x = year, y = catch * 1000, col = "reported", linetype = "reported"), linewidth = 1.2) +
  #geom_line(data = full_catch_df %>% filter(fmp_gear == "TRW") %>% group_by(year) %>% summarise(yearly_catch = sum(catch)), aes(x = year, y = yearly_catch, col = "reported", linetype = "reported"), linewidth = 1.2) +
  labs(x = "Year", y = "Catch", col = "", linetype = "") +
  theme_bw()

## save to data
saveRDS(fixed_gear_with_imputations_S1,file = file.path("Data", "5-Area", "fixed_gear_with_imputations_S1.RDS"))
saveRDS(fixed_gear_with_imputations_S2,file = file.path("Data", "5-Area", "fixed_gear_with_imputations_S2.RDS"))
saveRDS(trawl_gear_with_imputations_S1,file = file.path("Data", "5-Area", "trawl_gear_with_imputations_S1.RDS"))
saveRDS(trawl_gear_with_imputations_S2,file = file.path("Data", "5-Area", "trawl_gear_with_imputations_S2.RDS"))



