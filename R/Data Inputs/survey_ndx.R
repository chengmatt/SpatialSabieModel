#' 
#'  produce RPN's survey estimates
#'  for 3 Area and 5 Area model

library(tidyverse)
library(here)
library(RColorBrewer)
theme_set(theme_bw())

srv_data_df = read.csv(here("Data", "Survey", "Area RPNs for Strata 3 to 7.csv"))

## merge the eastern gulf regions
srv_data_df$NPFMC.Sablefish.Management.Area = ifelse(srv_data_df$NPFMC.Sablefish.Management.Area %in% c("West Yakutat","East Yakutat/Southeast"), "Eastern Gulf of Alaska", srv_data_df$NPFMC.Sablefish.Management.Area)

## like in the assessment drop US country prior to 1990
drop_ndx = srv_data_df$Country == 'United States' & srv_data_df$Year < 1990
srv_data_df = srv_data_df %>% filter(!drop_ndx)
## sum RPN and var_RPN to get area specific survey indicies
index_df = srv_data_df %>% group_by(Country, `NPFMC.Sablefish.Management.Area`, Year) %>% summarise(
  area_RPN = sum(RPN, na.rm = T), var_area_RPN = sum(RPN.Var, na.rm = T)
)
## only keep "exploitable areas"
index_exploitable_df = srv_data_df %>% filter(Exploitable == 1) %>% group_by(Country, `NPFMC.Sablefish.Management.Area`, Year) %>% summarise(
  area_RPN = sum(RPN, na.rm = T), var_area_RPN = sum(RPN.Var, na.rm = T)
)
index_df$type = "all"
index_exploitable_df$type = "exploit"
full_ndx = rbind(index_df, index_exploitable_df)
##
ggplot(index_df, aes(x = Year, y = area_RPN, col = Country)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~`NPFMC.Sablefish.Management.Area`, scales = "free_y") +
  theme_bw()
##
ggplot(index_exploitable_df, aes(x = Year, y = area_RPN, col = Country)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~`NPFMC.Sablefish.Management.Area`, scales = "free_y") +
  theme_bw()

## compare them on the same plot
ggplot(full_ndx, aes(x = Year, y = area_RPN, col = Country, linetype = type)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~`NPFMC.Sablefish.Management.Area`, scales = "free_y") +
  theme_bw()

## inspect the EGOA from US early, a weird drop off
EG = srv_data_df %>% filter(Country == "United States",`NPFMC.Sablefish.Management.Area` == "Eastern Gulf of Alaska")
ggplot(EG, aes(x = Year, y = RPN)) +
  geom_point(size = 1.1) +
  facet_wrap(~`Geographic.Area.Name`) +
  theme_bw() +
  ylim(0, NA)

EG %>% group_by(Year) %>% summarise(sum(RPN, na.rm = T))
EG %>% pivot_wider(id_cols = Year, names_from = `Geographic.Area.Name`, values_from = RPN)
## add estimates of uncertainty
## convert variance to a lognormal standard deviation
index_df$SE = log_sigma(sqrt(index_df$var_area_RPN) / index_df$area_RPN)
## calcualte Confidence intervals
CIs = lognormal_CI(expectation = index_df$area_RPN, sigma = index_df$SE, CI = 0.95)
index_df$LCI = CIs$lower
index_df$UCI = CIs$upper
##
ggplot(index_df, aes(x = Year, col = Country, fill = Country)) +
  geom_line(linewidth = 1.1, aes(y = area_RPN)) +
  geom_ribbon(aes(xmin = Year, xmax = Year, ymin = LCI, ymax = UCI), alpha = 0.4) +
  facet_wrap(~`NPFMC.Sablefish.Management.Area`, scales = "free_y") +
  theme_bw()

saveRDS(index_df, here("Data", "Survey", "regional_abundance_estimates.RDS"))

