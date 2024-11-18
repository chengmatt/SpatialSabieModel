# Purpose: To plot model diagnostic plots among 5-area models and 1-area models
# Creator: Matthew LH. Cheng
# Date: 10/17/24


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
library(compResidual)
source(here("R", 'Utility_Fxns.R'))

# Get Model Convergence Among 5-area models
out_5area_path = here("Output", "Final Models", "5-Area-1960")
out_1area_path = here("Output", "Final Models", "1-Area-1960")
colors = unname(ggthemes::ggthemes_data[["colorblind"]][["value"]]) # get colors

# Model Convergence -------------------------------------------------------
# Look at convergence real quick
files = c(list.files(out_5area_path), list.files(out_1area_path))
files = files[!str_detect(files, ".RDS")]

conv_all = data.frame()
for(i in 1:length(files)) {
  if(str_detect(files[i], "5-Area")) mod = readRDS(here(out_5area_path, files[i], 'sd_report.RDS'))
  if(str_detect(files[i], "1-Area")) mod = readRDS(here(out_1area_path, files[i], 'sd_report.RDS'))
  tmp = data.frame(mod = files[i], pdHess = mod$pdHess, grad = max(abs(mod$gradient.fixed)),
                   max_se = max(diag(mod$cov.fixed)))
  conv_all = rbind(tmp ,conv_all)
}


# Final 1-Area Comparisons + Diagnostics (No Tag vs. Tag) ------------------------------------------------------
model_list = list() # list to store models in 
model_sd = list() # list to store models in (sdrep)
model_dat = list() # list to store model data

# get model path f or 1 area models here
model_path = list(
  here(out_1area_path, "1-Area-1960-Final"),
  here(out_1area_path, "1-Area-1960-Final-Tagging")
)

model_name = c("1-Area", '1-Area-Tag') # set up model names

for(i in 1:length(model_path)) model_list[[i]] = readRDS(file.path(paste(model_path[[i]], "/mle_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_sd[[i]] = readRDS(file.path(paste(model_path[[i]], "/sd_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_dat[[i]] = readRDS(file.path(paste(model_path[[i]], "/data.RDS", sep = "")))

### Spawning Biomass Comparison --------------------------------------------------
ts_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get ssb
  ssb_tmp = model_sd[[i]]$value[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_se = model_sd[[i]]$sd[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_tmp_df = data.frame(year = 1960:2021, value = ssb_tmp, se = ssb_se, model = model_name[i], type = 'Sp Bio (kt)')
  # combine
  ts_all_df = rbind(ts_all_df, ssb_tmp_df)
} # end i

ssb_1area_plot = ggplot(ts_all_df, aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
                      color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  geom_ribbon(alpha = 0.2, color = NA) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = 'Model', fill = 'Model', lty = 'Model') 

ggsave(
  ssb_1area_plot,
  filename = here("figs", "Manuscript_Plots", "SSB_1Area_Comparison.png"),
  width = 10, height = 6
  )
  

### Model Diagnostics -------------------------------------------------------
##### Tag Recovery ------------------------------------------------------------
obs_tag_rec_1area = apply(model_list[[2]]$obs_tag_recovery, 3, sum)
pred_tag_rec_1area = apply(model_list[[2]]$pred_tag_recovery, 3, sum)
tag_rec_1area = data.frame(Obs = obs_tag_rec_1area, Pred = pred_tag_rec_1area, Year = 1:43, model = '1-Area-Tag') # combine for plotting

# Tag Recovery from 1-Area Model
ggplot(tag_rec_1area) +
  geom_point(aes(x = Year + 1978, y = Obs), size = 2) +
  geom_line(aes(x = Year + 1978, y = Pred), lwd = 1, alpha = 0.5) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = 'Model', fill = 'Model') 

##### Index -------------------------------------------------------------------
# Read in abundance index used
design_survey_index = readRDS(file = here("Data", "Survey", "regional_abundance_estimates.RDS"))
design_survey_index$area_lab = design_survey_index$NPFMC.Sablefish.Management.Area
design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "Aleutians" ~ "AI",
                                                                 area_lab == "Bering Sea" ~ "BS",
                                                                 area_lab == "Western Gulf of Alaska" ~ "WGOA",
                                                                 area_lab == "Central Gulf of Alaska" ~ "CGOA",
                                                                 area_lab == "Eastern Gulf of Alaska" ~ "EGOA",
                                                                 TRUE ~ area_lab)) %>% filter(Year != 2022)

# munge design based survey index
design_survey_index = design_survey_index %>% group_by(Country, Year) %>% 
  summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), 
            se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, 
            UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)

# Get predicted tag recovery
pred_srv_1area = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_tmp = reshape2::melt(model_list[[i]]$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, pred = value) %>%  
    mutate(model = model_name[i], nll = model_list[[i]]$nll[5])
  pred_srv_1area = rbind(pred_srv_1area, pred_srv_tmp)
} # end i

# residual munging
pred_srv_1area = pred_srv_1area %>% 
  mutate(Year = Year + 1959, 
         Country = ifelse(Country == 1, "Japan", "United States")) %>% 
  filter(Year >= 1979)

# Join datasets together
srv_bio_1area_df = design_survey_index %>%
  left_join(pred_srv_1area, by = c("Country", "Year")) %>% 
  mutate(resids = (log(sum_estimates) - log(pred)) / se,
         model = factor(model, levels = model_name))

# Pearson residuals
idx_resids = ggplot(srv_bio_1area_df) +
  geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
  geom_point(aes(x = Year, y = resids, color = model), size = 2, alpha = 0.75) +
  geom_smooth(aes(x = Year, y = resids, color = model, lty = model), se = FALSE, lwd = 1) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~Country, scales = "free") +
  labs(x = "Year", y = "Pearson Residuals", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.09,0.87),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

# Fits to index directly
idx_fits = ggplot(srv_bio_1area_df) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI)) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 1) +
  facet_wrap(~Country, scales = "free") +
  scale_color_manual(values = colors) +
  labs(x = "Year", y = "Relative Population Numbers", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.09,0.87),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

ggsave(
  plot_grid(idx_resids, idx_fits, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 20, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "Idx_Fits_1Area_Comparison.png"),
  width = 10, height = 10
)



###### Age and Length Comps ----------------------------------------------------
age_labels = c(
  paste("M", 2:31, sep = '-'),
  paste("F", 2:31, sep = '-')
) # sex-age labels

len_labels = c(
  paste("M", seq(41,99,2), sep = '-'),
  paste("F", seq(41,99,2), sep = '-')
) # sex-age labels

###### Fixed-Gear Age Comps ---------------------------------------------------------------
obs_fixed_age = reshape2::melt(model_list[[1]]$obs_fixed_catchatage) %>% 
  rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value) %>% 
  group_by(region, Year) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted fixed gear
pred_fixed_age = data.frame()
for(i in 1:length(model_list)) {
  pred_fixed_age_tmp = reshape2::melt(model_list[[i]]$pred_fixed_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[1])
  pred_fixed_age = rbind(pred_fixed_age, pred_fixed_age_tmp)
} # end i

fixed_age_1area_df = obs_fixed_age %>% 
  left_join(pred_fixed_age, by = c("Sex_Age", "Year")) 

# plot fixed gear age average fits
fixed_gear_avg = ggplot() +
  geom_col(fixed_age_1area_df %>% 
             filter(model == model_name[1]) %>% 
             group_by(Sex_Age, model) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(fixed_age_1area_df %>% 
            group_by(Sex_Age, model) %>%
            summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = age_labels[seq(1,60, 4)], breaks = seq(1, 60, 4)) +
  scale_color_manual(values = colors) +
  labs(x = "Sex-Age Category", y = "Fixed-Gear Age Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) + 
  theme(legend.position = c(0.93,0.86),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

fixed_gear_avg

# plot residuals 
# Get OSA residuals first
osa_all = data.frame()
for(i in 1:length(model_list)) { 
  obs = t(model_list[[i]]$obs_fixed_catchatage[,1,40:62]) # get observed
  pred = t(model_list[[i]]$pred_fixed_catchatage[,1,40:62]) # get predicted
  rownames(pred) = 1999:2021 # some row name munging
  osa_res_tmp = get_osa_res(obs = obs, pred = pred, iss = 1, iter_wt = 1, index = age_labels, drop_bin = 2) %>% 
    mutate(model = model_name[i]) # get OSA residuals
  osa_all = rbind(osa_all, osa_res_tmp)
} # end i

# OSA residuals
osa_1area = ggplot(osa_all, aes(Year, as.numeric(Sex_Age), size=abs(resid), color=resid>0)) + 
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Age Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_wrap(~model) +
  scale_y_continuous(labels = age_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

osa_1area

ggsave(
  plot_grid(osa_1area, fixed_gear_avg, ncol = 1,
            align = "hv", axis = 'l', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "FixedAge_1Area_Comparison.png"),
  width = 10, height = 10
)

 
###### Survey Age Comps --------------------------------------------------------
obs_srv_age = reshape2::melt(model_list[[1]]$obs_srv_catchatage) %>% 
  rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value, Country = Var4) %>% 
  group_by(region, Year, Country) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted srv age
pred_srv_age = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_age_tmp = reshape2::melt(model_list[[i]]$pred_srv_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Country = Var4, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[4])
  pred_srv_age = rbind(pred_srv_age, pred_srv_age_tmp)
} # end i

srv_age_1area_df = obs_srv_age %>% 
  left_join(pred_srv_age, by = c("Sex_Age", "region", "Year", "Country")) %>% 
  mutate(
    Country = ifelse(Country == 1, "Japan", "United States"),
    model = factor(model, levels = model_name)
  ) 

# plot fixed gear age average fits
srv_gear_avg = ggplot() +
  geom_col(srv_age_1area_df %>%
             filter(model == model_name[1]) %>%
             group_by(Sex_Age, model, Country) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(srv_age_1area_df %>%
              group_by(Sex_Age, model, Country) %>%
              summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = age_labels[seq(1,60, 4)], breaks = seq(1, 60, 4)) +
  scale_color_manual(values = colors) +
  facet_wrap(~Country) +
  labs(x = "Sex-Age Category", y = "Survey-Gear Age Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9,0.86),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))

srv_gear_avg

# plot residuals
# Get OSA residuals first
osa_all = data.frame()
for(i in 1:length(model_list)) {
  obs_jp = t(model_list[[i]]$obs_srv_catchatage[,1,c(22, seq(26, 34, 2)),1]) # get observed
  pred_jp = t(model_list[[i]]$pred_srv_catchatage[,1,c(22, seq(26, 34, 2)),1]) # get predicted
  obs_dom = t(model_list[[i]]$obs_srv_catchatage[,1,37:62,2]) # get observed
  pred_dom = t(model_list[[i]]$pred_srv_catchatage[,1,37:62,2]) # get predicted
  rownames(pred_jp) = c(1981, seq(1985, 1993, 2)) # some row name munging
  rownames(pred_dom) = 1996:2021 # some row name munging
  # get OSA Japan
  osa_res_jp_tmp = get_osa_res(obs = obs_jp, pred = pred_jp, iss = 1, iter_wt = 1, 
                               index = age_labels, drop_bin = 1) %>% mutate(model = model_name[i], Country = 'Japan') 
  # get OSA US
  osa_res_us_tmp = get_osa_res(obs = obs_dom, pred = pred_dom, iss = 1, iter_wt = 1, 
                               index = age_labels, drop_bin = 1) %>%  mutate(model = model_name[i], Country = 'United States') 
  osa_all = rbind(osa_all, osa_res_jp_tmp, osa_res_us_tmp)
} # end i

# OSA residuals
osa_1area = ggplot(osa_all, aes(Year, as.numeric(Sex_Age), size=abs(resid), color=resid>0)) +
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Age Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_grid(model~Country) +
  scale_y_continuous(labels = age_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  plot_grid(osa_1area, srv_gear_avg, ncol = 1,
            align = "hv", axis = 'l', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "SrvAge_1Area_Comparison.png"),
  width = 15, height = 12
)


###### Trawl Length Comps --------------------------------------------------------
obs_trawl_len = reshape2::melt(model_list[[1]]$obs_trwl_catchatlgth) %>% 
  rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
  group_by(region, Year) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted fixed gear
pred_trawl_len = data.frame()
for(i in 1:length(model_list)) {
  pred_trawl_len_tmp = reshape2::melt(model_list[[i]]$pred_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[1])
  pred_trawl_len = rbind(pred_trawl_len, pred_trawl_len_tmp)
} # end i

trawl_len_1area_df = obs_trawl_len %>% 
  left_join(pred_trawl_len, by = c("Sex_Len", "Year")) 

# plot fixed gear age average fits
trawl_gear_avg = ggplot() +
  geom_col(trawl_len_1area_df %>% 
             filter(model == model_name[1]) %>% 
             group_by(Sex_Len, model) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(trawl_len_1area_df %>% 
              group_by(Sex_Len, model) %>%
              summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = len_labels[seq(1,60, 4)], breaks = seq(1, 60, 4)) +
  scale_color_manual(values = colors) +
  labs(x = "Sex-Length Category", y = "Trawl-Gear Length Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) + 
  theme(legend.position = c(0.93,0.86),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

trawl_gear_avg

# plot residuals 
# Get OSA residuals first
osa_all = data.frame()
for(i in 1:length(model_list)) { 
  obs = t(model_list[[i]]$obs_trwl_catchatlgth[,1,40:62]) # get observed
  pred = t(model_list[[i]]$pred_trwl_catchatlgth[,1,40:62]) # get predicted
  rownames(pred) = 1999:2021 # some row name munging
  osa_res_tmp = get_osa_res(obs = obs, pred = pred, iss = 1, iter_wt = 1, index = age_labels, drop_bin = 1) %>% 
    mutate(model = model_name[i]) %>% rename(Sex_Len = Sex_Age) # get OSA residuals
  osa_all = rbind(osa_all, osa_res_tmp)
} # end i

# OSA residuals
osa_1area = ggplot(osa_all, aes(Year, as.numeric(Sex_Len), size=abs(resid), color=resid>0)) + 
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Length Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_wrap(~model) +
  scale_y_continuous(labels = len_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

osa_1area

ggsave(
  plot_grid(osa_1area, trawl_gear_avg, ncol = 1,
            align = "hv", axis = 'l', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "TrwlLen_1Area_Comparison.png"),
  width = 10, height = 10
)


# Final 5-Area Model Diagnostics ------------------------------------------------------
model_list = list() # list to store models in 
model_sd = list() # list to store models in (sdrep)
model_dat = list() # list to store model data

# get model path f or 1 area models here
model_path = list(
  here(out_5area_path, "5-Area-1960-03-FishBlock")
)

model_name = c("5-Area") # set up model names

for(i in 1:length(model_path)) model_list[[i]] = readRDS(file.path(paste(model_path[[i]], "/mle_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_sd[[i]] = readRDS(file.path(paste(model_path[[i]], "/sd_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_dat[[i]] = readRDS(file.path(paste(model_path[[i]], "/data.RDS", sep = "")))

### Spawning Biomass Comparison --------------------------------------------------
ts_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get ssb
  ssb_tmp = model_sd[[i]]$value[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_se = model_sd[[i]]$sd[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_tmp_df = data.frame(year = 1960:2021, value = ssb_tmp, se = ssb_se, model = model_name[i], type = 'Sp Bio (kt)',
                          region = rep(c("BS","AI","WGOA","CGOA","EGOA"), each = length(1960:2021))) %>% 
    mutate(region = factor(region, levels = c("BS","AI","WGOA","CGOA","EGOA")))
  # combine
  ts_all_df = rbind(ts_all_df, ssb_tmp_df)
} # end i

ssb_5area_plot = ggplot(ts_all_df, aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
                      color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_wrap(~region, nrow = 1) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = 'Model', fill = 'Model', lty = 'Model') 

ggsave(
  ssb_5area_plot,
  filename = here("figs", "Manuscript_Plots", "SSB_5Area_Comparison.png"),
  width = 14, height = 6
)

### Model Diagnostics -------------------------------------------------------
##### Tag Recovery ------------------------------------------------------------
obs_tag_rec = reshape2::melt(model_list[[1]]$obs_tag_recovery) %>% 
  rename(Age = Var1, Event = Var2, Region = Var3, Year = Var4, Obs = value)
pred_tag_rec = reshape2::melt(model_list[[1]]$pred_tag_recovery) %>% 
  rename(Age = Var1, Event = Var2, Region = Var3, Year = Var4, Pred = value)

tag_rec = cbind(obs_tag_rec, Pred = pred_tag_rec$Pred) %>% 
  mutate(
    Region = case_when(
      Region == 1 ~ "BS",
      Region == 2 ~ "AI",
      Region == 3 ~ "WGOA",
      Region == 4 ~ "CGOA",
      Region == 5 ~ "EGOA"
    ), Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# Tag Recovery from 5-Area Model
area_5_tagrec_plot = ggplot(tag_rec %>% group_by(Year, Region) %>% summarize(Obs = sum(Obs), Pred = sum(Pred))) +
  geom_point(mapping = aes(x = Year + 1978, y = Obs), size = 2) +
  geom_line(mapping = aes(x = Year + 1978, y = Pred), lwd = 1, alpha = 0.8) +
  facet_wrap(~Region, nrow = 1) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = 'Model', fill = 'Model') 

ggsave(
  area_5_tagrec_plot,
  filename = here("figs", "Manuscript_Plots", "TagRecFits_Comp5area.png"),
  width = 15, height = 8
)

# Tag Recovery Comparison w/ 1-Area Tag Model
tag_area_comparison = rbind(
  tag_rec %>% group_by(Year) %>% summarize(Obs = sum(Obs), Pred = sum(Pred)) %>% mutate(model = '5-Area'),
  tag_rec_1area
)

# Tag Recovery from 1-Area Model comparison
tag_rec_comp_plot = ggplot() +
  geom_point(tag_area_comparison %>% filter(model == '5-Area'), mapping = aes(x = Year + 1978, y = Obs), size = 2) +
  geom_line(tag_area_comparison, mapping = aes(x = Year + 1978, y = Pred, color = model, lty = model), lwd = 1, alpha = 1) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = 'Model', fill = 'Model', lty = 'Model') 

ggsave(
  tag_rec_comp_plot,
  filename = here("figs", "Manuscript_Plots", "TagRecFits_Comp5A_1A.png"),
  width = 10, height = 10
)

##### Index -------------------------------------------------------------------
# Read in abundance index used
design_survey_index = readRDS(file = here("Data", "Survey", "regional_abundance_estimates.RDS"))
design_survey_index$area_lab = design_survey_index$NPFMC.Sablefish.Management.Area
design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "Aleutians" ~ "AI",
                                                                 area_lab == "Bering Sea" ~ "BS",
                                                                 area_lab == "Western Gulf of Alaska" ~ "WGOA",
                                                                 area_lab == "Central Gulf of Alaska" ~ "CGOA",
                                                                 area_lab == "Eastern Gulf of Alaska" ~ "EGOA",
                                                                 TRUE ~ area_lab)) %>% filter(Year != 2022)

# munge design based survey index
design_survey_index = design_survey_index %>% group_by(area_lab, Country, Year) %>% 
  summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), 
            se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, 
            UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)

# Get predicted idx
pred_srv = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_tmp = reshape2::melt(model_list[[i]]$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, pred = value) %>%  
    mutate(model = model_name[i], nll = model_list[[i]]$nll[5])
  pred_srv = rbind(pred_srv, pred_srv_tmp)
} # end i

# residual munging
pred_srv = pred_srv %>% 
  mutate(Year = Year + 1959, 
         area_lab = case_when(
           area_lab == 1 ~ "BS",
           area_lab == 2 ~ "AI",
           area_lab == 3 ~ "WGOA",
           area_lab == 4 ~ "CGOA",
           area_lab == 5 ~ "EGOA"),
         Country = ifelse(Country == 1, "Japan", "United States")) %>% 
  filter(Year >= 1979)

# Join datasets together
srv_bio_df = design_survey_index %>%
  left_join(pred_srv, by = c("Country", "area_lab", "Year")) %>% 
  mutate(resids = (log(sum_estimates) - log(pred)) / se,
         model = factor(model, levels = model_name),
         area_lab = factor(area_lab, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# Pearson residuals
idx_resids = ggplot(srv_bio_df) +
  geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
  geom_point(aes(x = Year, y = resids, color = model), size = 2, alpha = 0.75) +
  geom_smooth(aes(x = Year, y = resids, color = model, lty = model), se = FALSE, lwd = 1) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_grid(area_lab~Country, scales = "free") +
  labs(x = "Year", y = "Pearson Residuals", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

# Fits to index directly
idx_fits = ggplot(srv_bio_df) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 1) +
  facet_wrap(~Country, scales = "free") +
  scale_color_manual(values = colors) +
  facet_grid(area_lab~Country, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

# plot aggregated fits and compare to 1 area
srv_5area_agg_preds = srv_bio_df %>% 
  group_by(Country,Year) %>% 
  summarise(sum_estimates = sum(sum_estimates, na.rm = T), sum_var = sum(sum_var, na.rm = T), 
            se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, 
            UCI = lognormal_CI(sum_estimates, se, 0.95)$upper, pred = sum(pred)) %>% 
  mutate(model = '5-Area')

# Combine 1-area and 5-area predictions
srv_bio_agg_preds = srv_bio_1area_df %>% 
  filter(model == '1-Area') %>% 
  select(Country, Year, sum_estimates, sum_var, se, LCI, UCI, pred, model) %>% 
  rbind(srv_5area_agg_preds) %>% 
  mutate(resids = (log(sum_estimates) - log(pred)) / se)

# Combined aggregated fits comparing 5-area and 1-area model
comb_idx_fits = ggplot() +
  geom_pointrange(srv_bio_agg_preds, mapping = aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI)) +
  geom_line(srv_bio_agg_preds, mapping = aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 1) +
  facet_wrap(~Country, scales = 'free') +
  scale_color_manual(values = colors) +
  labs(x = "Year", y = "Relative Population Numbers", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.6, 0.83),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

comb_idx_fits

# combined residuals
comb_idx_resids = ggplot(srv_bio_agg_preds) +
  geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
  geom_point(aes(x = Year, y = resids, color = model), size = 2, alpha = 0.75) +
  geom_smooth(aes(x = Year, y = resids, color = model, lty = model), se = FALSE, lwd = 1) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~Country, scales = "free") +
  labs(x = "Year", y = "Pearson Residuals", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

comb_idx_resids

# plot 5-area region specific comparisons
ggsave(
  plot_grid(idx_resids, idx_fits, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            align = 'v', axis = 'l',
            label_size = 20, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "Idx_Fits_5Area_Comparison.png"),
  width = 10, height = 13
)

# plot combined comparisons
ggsave(
  plot_grid(comb_idx_resids, comb_idx_fits, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            align = 'v', axis = 'l',
            label_size = 20, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "Idx_Fits_5A_1A_Comp.png"),
  width = 13, height = 10
)

###### Age and Length Comps ----------------------------------------------------
###### Fixed-Gear Age Comps ----------------------------------------------------
obs_fixed_age = reshape2::melt(model_list[[1]]$obs_fixed_catchatage) %>% 
  rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value) %>% 
  group_by(region, Year) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted fixed gear
pred_fixed_age = data.frame()
for(i in 1:length(model_list)) {
  pred_fixed_age_tmp = reshape2::melt(model_list[[i]]$pred_fixed_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[1])
  pred_fixed_age = rbind(pred_fixed_age, pred_fixed_age_tmp)
} # end i

fixed_age_df = obs_fixed_age %>% 
  left_join(pred_fixed_age, by = c("Sex_Age", "Year", 'region')) %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# plot fixed gear age average fits
fixed_gear_avg = ggplot() +
  geom_col(fixed_age_df %>% 
             filter(model == model_name[1]) %>% 
             group_by(Sex_Age, model, region) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(fixed_age_df %>% 
              group_by(Sex_Age, model, region) %>%
              summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = age_labels[seq(1,60, 5)], breaks = seq(1, 60, 5)) +
  scale_color_manual(values = colors) +
  facet_wrap(~region, nrow = 1) +
  labs(x = "Sex-Age Category", y = "Fixed-Gear Age Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) + 
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

fixed_gear_avg

# plot residuals 
# Get OSA residuals first
n_regions = 5
osa_all = data.frame()
for(i in 1:length(model_list)) { 
  for(r in 1:n_regions) {
    obs = t(model_list[[i]]$obs_fixed_catchatage[,r,40:62]) # get observed
    pred = t(model_list[[i]]$pred_fixed_catchatage[,r,40:62]) # get predicted
    rownames(pred) = 1999:2021 # some row name munging
    osa_res_tmp = get_osa_res(obs = obs, pred = pred, iss = 1, iter_wt = 1, index = age_labels, drop_bin = 1) %>% 
      mutate(model = model_name[i], Region = r) # get OSA residuals
    osa_all = rbind(osa_all, osa_res_tmp)
    
  } # end r
} # end i

# name regions
osa_all = osa_all %>% 
  mutate(
    Region = case_when(
      Region == 1 ~ "BS",
      Region == 2 ~ "AI",
      Region == 3 ~ "WGOA",
      Region == 4 ~ "CGOA",
      Region == 5 ~ "EGOA"
    ), Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# OSA residuals
osa_5area = ggplot(osa_all, aes(Year, as.numeric(Sex_Age), size=abs(resid), color=resid>0)) + 
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Age Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_wrap(~Region, nrow = 1) +
  scale_y_continuous(labels = age_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

osa_5area

ggsave(
  plot_grid(osa_5area, fixed_gear_avg, ncol = 1,
            align = "v", axis = 'lr', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1),
  filename = here("figs", "Manuscript_Plots", "FixedAge_5Area_Comparison.png"),
  width = 15, height = 13
)

###### Survey Age Comps --------------------------------------------------------
obs_srv_age = reshape2::melt(model_list[[1]]$obs_srv_catchatage) %>% 
  rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value, Country = Var4) %>% 
  group_by(region, Year, Country) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted srv age
pred_srv_age = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_age_tmp = reshape2::melt(model_list[[i]]$pred_srv_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Country = Var4, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[4])
  pred_srv_age = rbind(pred_srv_age, pred_srv_age_tmp)
} # end i

srv_age_df = obs_srv_age %>% 
  left_join(pred_srv_age, by = c("Sex_Age", "region", "Year", "Country")) %>% 
  mutate(
    Country = ifelse(Country == 1, "Japan", "United States"),
    model = factor(model, levels = model_name),
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# plot fixed gear age average fits
srv_gear_avg = ggplot() +
  geom_col(srv_age_df %>%
             filter(model == model_name[1]) %>%
             group_by(Sex_Age, model, Country, region) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(srv_age_df %>%
              group_by(Sex_Age, model, Country, region) %>%
              summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = age_labels[seq(1,60, 4)], breaks = seq(1, 60, 4)) +
  scale_color_manual(values = colors) +
  facet_grid(Country~region) +
  labs(x = "Sex-Age Category", y = "Survey-Gear Age Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.background = element_rect(fill = "transparent", colour = NA))


# plot residuals
# Get OSA residuals first
osa_all = data.frame()
for(i in 1:length(model_list)) {
  for(r in 1:n_regions) {
    obs_jp = t(model_list[[i]]$obs_srv_catchatage[,r,c(22, seq(26, 34, 2)),1]) # get observed
    pred_jp = t(model_list[[i]]$pred_srv_catchatage[,r,c(22, seq(26, 34, 2)),1]) # get predicted
    obs_dom = t(model_list[[i]]$obs_srv_catchatage[,r,37:62,2]) # get observed
    pred_dom = t(model_list[[i]]$pred_srv_catchatage[,r,37:62,2]) # get predicted
    rownames(pred_jp) = c(1981, seq(1985, 1993, 2)) # some row name munging
    
    if(sum(which(rowSums(obs_dom) == 0)) > 0) {
      pred_dom = pred_dom[-which(rowSums(obs_dom) == 0), ]
      rownames(pred_dom) = (1996:2021)[-which(rowSums(obs_dom) == 0)] # some row name munging
      obs_dom = obs_dom[-which(rowSums(obs_dom) == 0), ]
    } else {
      rownames(pred_dom) = (1996:2021)
    } # if else for 0s in rowSums

    # get OSA Japan
    osa_res_jp_tmp = get_osa_res(obs = obs_jp, pred = pred_jp, iss = 1, iter_wt = 1, 
                                 index = age_labels, drop_bin = 1) %>% mutate(model = model_name[i], Country = 'Japan', Region = r) 
    # get OSA US
    osa_res_us_tmp = get_osa_res(obs = obs_dom, pred = pred_dom, iss = 1, iter_wt = 1, 
                                 index = age_labels, drop_bin = 1) %>%  mutate(model = model_name[i], Country = 'United States', Region = r) 
    osa_all = rbind(osa_all, osa_res_jp_tmp, osa_res_us_tmp)
  } # end r
} # end i

# name regions
osa_all = osa_all %>% 
  mutate(
    Region = case_when(
      Region == 1 ~ "BS",
      Region == 2 ~ "AI",
      Region == 3 ~ "WGOA",
      Region == 4 ~ "CGOA",
      Region == 5 ~ "EGOA"
    ), Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# OSA residuals
osa_5area = ggplot(osa_all, aes(Year, as.numeric(Sex_Age), size=abs(resid), color=resid>0)) +
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Age Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_grid(Country~Region) +
  scale_y_continuous(labels = age_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))

osa_5area

ggsave(
  plot_grid(osa_5area, srv_gear_avg, ncol = 1,
            align = "v", axis = 'lr', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1),
  filename = here("figs", "Manuscript_Plots", "SrvAge_5Area_Comparison.png"),
  width = 17, height = 13
)


###### Trawl Length Comps --------------------------------------------------------
obs_trawl_len = reshape2::melt(model_list[[1]]$obs_trwl_catchatlgth) %>% 
  rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
  group_by(region, Year) %>% 
  mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
  filter(total != 0)

# get predicted fixed gear
pred_trawl_len = data.frame()
for(i in 1:length(model_list)) {
  pred_trawl_len_tmp = reshape2::melt(model_list[[i]]$pred_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[1])
  pred_trawl_len = rbind(pred_trawl_len, pred_trawl_len_tmp)
} # end i

trawl_len_df = obs_trawl_len %>% 
  left_join(pred_trawl_len, by = c("Sex_Len", "Year", 'region')) %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# plot fixed gear age average fits
trawl_gear_avg = ggplot() +
  geom_col(trawl_len_df %>% 
             filter(model == model_name[1]) %>% 
             group_by(Sex_Len, model, region) %>%
             summarize(Obs_Mean = mean(obs_prop)),
           mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.25, color = 'black') + # observed data
  geom_line(trawl_len_df %>% 
              group_by(Sex_Len, model, region) %>%
              summarize(Pred_Mean = mean(pred)), # predictions
            mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
  scale_x_continuous(labels = len_labels[seq(1,60, 4)], breaks = seq(1, 60, 4)) +
  scale_color_manual(values = colors) +
  facet_wrap(~region, nrow = 1) +
  labs(x = "Sex-Length Category", y = "Trawl-Gear Length Compositions", color = "Model", lty = "Model") +
  theme_bw(base_size = 15) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

trawl_gear_avg

# plot residuals 
# Get OSA residuals first
osa_all = data.frame()
for(i in 1:length(model_list)) { 
  for(r in 1:n_regions) {
    obs = t(model_list[[i]]$obs_trwl_catchatlgth[,r,40:62]) # get observed
    pred = t(model_list[[i]]$pred_trwl_catchatlgth[,r,40:62]) # get predicted
    rownames(pred) = 1999:2021 # some row name munging
    osa_res_tmp = get_osa_res(obs = obs, pred = pred, iss = 1, iter_wt = 1, index = age_labels, drop_bin = 1) %>% 
      mutate(model = model_name[i], Region = r) %>% rename(Sex_Len = Sex_Age) # get OSA residuals
    osa_all = rbind(osa_all, osa_res_tmp)
  } # end r
} # end i

# name regions
osa_all = osa_all %>% 
  mutate(
    Region = case_when(
      Region == 1 ~ "BS",
      Region == 2 ~ "AI",
      Region == 3 ~ "WGOA",
      Region == 4 ~ "CGOA",
      Region == 5 ~ "EGOA"
    ), Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))  

# OSA residuals
osa_5area = ggplot(osa_all, aes(Year, as.numeric(Sex_Len), size=abs(resid), color=resid>0)) + 
  geom_point(alpha = 0.35) +
  labs(x = "Year", y = 'Sex-Length Category', size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
  guides(color = guide_legend(order = 0), size = guide_legend(order = 1)) +
  facet_wrap(~Region, nrow = 1) +
  scale_y_continuous(labels = len_labels[-1][seq(2,60, 4)], breaks = seq(2, 60, 4)) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "top",
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) 

osa_5area

ggsave(
  plot_grid(osa_5area, trawl_gear_avg, ncol = 1,
            align = "v", axis = 'l', rel_heights = c(0.6, 0.4),
            labels = c('A', 'B'),
            label_size = 25, hjust = -1.85),
  filename = here("figs", "Manuscript_Plots", "TrwlLen_5Area_Comparison.png"),
  width = 15, height = 13
)


# 5-Area Model Comparisons ------------------------------------------------
colors = ggthemes::tableau_color_pal('Classic 10 Medium')(10)

# Filter to converged models
models_5area = conv_all %>% 
  filter(grad < 0.005, str_detect(mod, '5-Area'), pdHess == TRUE,
         !is.nan(max_se), max_se < 100)

# Extract models
model_list = list() # list to store models in 
model_sd = list() # list to store models in (sdrep)
model_dat = list() # list to store model data
model_path = list()

# get model path for 5 area models here
for(i in 1:nrow(models_5area)) model_path[[i]] = here(out_5area_path, models_5area$mod[i])
for(i in 1:length(model_path)) model_list[[i]] = readRDS(file.path(paste(model_path[[i]], "/mle_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_sd[[i]] = readRDS(file.path(paste(model_path[[i]], "/sd_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_dat[[i]] = readRDS(file.path(paste(model_path[[i]], "/data.RDS", sep = "")))

# Set up model names
model_name = str_remove(models_5area$mod, "5-Area-1960-")
model_name = sub("_", "",sub("_", "", model_name)) # remove underscores...

# for factoring
model_levels = c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move",
                 "02-NegBin", "03-Decadal", "03-FishBlock",
                 "03-Space", "04-SptQ")

### Spawning Biomass Comparison --------------------------------------------------
ts_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get ssb
  ssb_tmp = model_sd[[i]]$value[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_se = model_sd[[i]]$sd[names(model_sd[[i]]$value) == 'SSB_yr']
  ssb_tmp_df = data.frame(year = 1960:2021, value = ssb_tmp, se = ssb_se,
                          cv = ssb_se / ssb_tmp, model = model_name[i], type = 'Sp Bio (kt)',
                          region = rep(c("BS","AI","WGOA","CGOA","EGOA"), each = length(1960:2021))) %>% 
    mutate(region = factor(region, levels = c("BS","AI","WGOA","CGOA","EGOA")))
  # combine
  ts_all_df = rbind(ts_all_df, ssb_tmp_df)
} # end i

# relevel model names
ts_all_df = ts_all_df %>%
  mutate(model = factor(model, levels = model_levels),
         # Creating facets to combine comparisons
         movement = ifelse(model %in% c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
         taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
         tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
         sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA)
           )

# Plot SSB movement comparisons
move_ssb = ggplot(ts_all_df %>% filter(!is.na(movement)), 
       aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
           color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_grid(movement~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB tag likelihood comparisons
taglike_ssb = ggplot(ts_all_df %>% filter(!is.na(taglike)), 
       aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
           color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  scale_fill_manual(values = colors[-c(1,2,3)]) +
  facet_grid(taglike~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB tag reporting rates 
tagrep_ssb = ggplot(ts_all_df %>% filter(!is.na(tagrep)), 
       aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
           color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB spatial Q
sptq_ssb = ggplot(ts_all_df %>% filter(!is.na(sptq)), 
       aes(x = year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se, 
           color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Sp Bio (kt)', color = '', fill = '', lty = '') 


ggsave(
  plot_grid(move_ssb, taglike_ssb, tagrep_ssb, sptq_ssb, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18),
  filename = here("figs", "Manuscript_Plots", "SSB_5Area_Comparison.png"),
  width = 15, height = 10
) 


### Spawning Biomass Comparison (CV) -----------------------------------------------------
# Plot SSB movement comparisons
move_ssb = ggplot(ts_all_df %>% filter(!is.na(movement)), 
                  aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_grid(movement~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'CV Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB tag likelihood comparisons
taglike_ssb = ggplot(ts_all_df %>% filter(!is.na(taglike)), 
                     aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  scale_fill_manual(values = colors[-c(1,2,3)]) +
  facet_grid(taglike~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'CV Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB tag reporting rates 
tagrep_ssb = ggplot(ts_all_df %>% filter(!is.na(tagrep)), 
                    aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'CV Sp Bio (kt)', color = '', fill = '', lty = '') 

# Plot SSB spatial Q
sptq_ssb = ggplot(ts_all_df %>% filter(!is.na(sptq)), 
                  aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'CV Sp Bio (kt)', color = '', fill = '', lty = '') 

ggsave(
  plot_grid(move_ssb, taglike_ssb, tagrep_ssb, sptq_ssb, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18),
  filename = here("figs", "Manuscript_Plots", "SSB_CV_5Area_Comparison.png"),
  width = 15, height = 10
) 

### Recruitment Comparison --------------------------------------------------
ts_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get ssb
  rec_tmp = model_sd[[i]]$value[names(model_sd[[i]]$value) == 'recruitment_yr']
  rec_se = model_sd[[i]]$sd[names(model_sd[[i]]$value) == 'recruitment_yr']
  rec_tmp_df = data.frame(year = 1960:2021, value = rec_tmp, se = rec_se, 
                          model = model_name[i], type = 'Recrmt', cv = rec_se / rec_tmp, 
                          region = rep(c("BS","AI","WGOA","CGOA","EGOA"), each = length(1960:2021))) %>% 
    mutate(region = factor(region, levels = c("BS","AI","WGOA","CGOA","EGOA")))
  # combine
  ts_all_df = rbind(ts_all_df, rec_tmp_df)
} # end i

# relevel model names
ts_all_df = ts_all_df %>%
  mutate(model = factor(model, levels = model_levels),
         # Creating facets to combine comparisons
         movement = ifelse(model %in% c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
         taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
         tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
         sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA)
  )

# Plot rec movement comparisons
move_rec = ggplot(ts_all_df %>% filter(!is.na(movement)), 
                  aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_grid(movement~region) +
  coord_cartesian(ylim = c(0,160)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec tag likelihood comparisons
taglike_rec = ggplot(ts_all_df %>% filter(!is.na(taglike)), 
                     aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  scale_fill_manual(values = colors[-c(1,2,3)]) +
  facet_grid(taglike~region) +
  coord_cartesian(ylim = c(0,160)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec tag reporting rates 
tagrep_rec = ggplot(ts_all_df %>% filter(!is.na(tagrep)), 
                    aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~region) +
  coord_cartesian(ylim = c(0,160)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec spatial Q
sptq_rec = ggplot(ts_all_df %>% filter(!is.na(sptq)), 
                  aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~region) +
  coord_cartesian(ylim = c(0,160)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Recrmt (millions)", color = '', fill = '', lty = '') 


ggsave(
  plot_grid(move_rec, taglike_rec, tagrep_rec, sptq_rec, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18),
  filename = here("figs", "Manuscript_Plots", "Rec_5Area_Comparison.png"),
  width = 15, height = 10
) 


### Recruitment Comparison (CV) ---------------------------------------------
# Plot rec movement comparisons
move_rec = ggplot(ts_all_df %>% filter(!is.na(movement)), 
                  aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_grid(movement~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.135, 0.78),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "CV Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec tag likelihood comparisons
taglike_rec = ggplot(ts_all_df %>% filter(!is.na(taglike)), 
                     aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  scale_fill_manual(values = colors[-c(1,2,3)]) +
  facet_grid(taglike~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.135, 0.78),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "CV Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec tag reporting rates 
tagrep_rec = ggplot(ts_all_df %>% filter(!is.na(tagrep)), 
                    aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.135, 0.78),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "CV Recrmt (millions)", color = '', fill = '', lty = '') 

# Plot rec spatial Q
sptq_rec = ggplot(ts_all_df %>% filter(!is.na(sptq)), 
                  aes(x = year, y = cv, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.135, 0.78),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "CV Recrmt (millions)", color = '', fill = '', lty = '') 


ggsave(
  plot_grid(move_rec, taglike_rec, tagrep_rec, sptq_rec, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18, label_y = 1.035),
  filename = here("figs", "Manuscript_Plots", "Rec_CV_5Area_Comparison.png"),
  width = 15, height = 10
) 


### Mean Female Age ---------------------------------------------------------
ts_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get mean female age
  tmp = reshape2::melt(model_list[[i]]$natage_f) %>% mutate(type = 'Mean Female Age') %>% 
    group_by(Var3, Var2, type) %>% # group by region and year
    mutate(prop = value / sum(value)) %>% 
    summarize(value = sum(prop * 2:31)) %>% 
    filter(Var3 != 63) %>% 
    rename(Var1 = Var3) %>% 
    rename(year = Var1, region = Var2) %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), 
      region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
      year = year + 1959,
      model = model_name[i]) 
  
  # combine
  ts_all_df = rbind(ts_all_df, tmp)
} # end i

# relevel model names
ts_all_df = ts_all_df %>%
  mutate(model = factor(model, levels = model_levels),
         # Creating facets to combine comparisons
         movement = ifelse(model %in% c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
         taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
         tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
         sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA)
  )

# Plot mean female age movement comparisons
move_age = ggplot(ts_all_df %>% filter(!is.na(movement)), 
                  aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_grid(movement~region) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.29),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Mean Female Age", color = '', fill = '', lty = '') 

# Plot female age tag likelihood comparisons
taglike_age = ggplot(ts_all_df %>% filter(!is.na(taglike)), 
                     aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  scale_fill_manual(values = colors[-c(1,2,3)]) +
  facet_grid(taglike~region) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.29),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Mean Female Age", color = '', fill = '', lty = '') 

# Plot female age tag reporting rates 
tagrep_age = ggplot(ts_all_df %>% filter(!is.na(tagrep)), 
                    aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~region) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.29),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Mean Female Age", color = '', fill = '', lty = '') 

# Plot female age spatial Q
sptq_age = ggplot(ts_all_df %>% filter(!is.na(sptq)), 
                  aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  scale_fill_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~region) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.95, 0.29),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = "Mean Female Age", color = '', fill = '', lty = '') 


ggsave(
  plot_grid(move_age, taglike_age, tagrep_age, sptq_age, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18),
  filename = here("figs", "Manuscript_Plots", "MeanFemaleAge_5Area_Comparison.png"),
  width = 15, height = 13
) 

### Movement  ---------------------------------------------------------------
move_all_df = data.frame()
for(i in 1:length(model_path)) {
  # get movement information
  move_tmp = reshape2::melt(model_list[[i]]$movement_matrix)
  move_se = model_sd[[i]]$sd[names(model_sd[[i]]$value) == 'movement_matrix']
  
  # create se, cv and model name for dataframe
  move_tmp$se = move_se
  move_tmp$cv = move_tmp$se / move_tmp$value
  move_tmp$model = model_name[i]
  
  # combine
  move_all_df = rbind(move_all_df, move_tmp)
} # end i

# Do some quick munging
move_all_df = move_all_df %>% 
  rename(From = Var1, To = Var2, Time = Var3, Age = Var4) %>% 
  mutate(
    From = case_when(
      From == 1 ~ "BS",
      From == 2 ~ "AI",
      From == 3 ~ "WGOA",
      From == 4 ~ "CGOA",
      From == 5 ~ "EGOA"
    ), 
    From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), 
    To = factor(To, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    model = factor(model, levels = model_levels),
    # Creating facets to combine comparisons
    movement = ifelse(model %in% c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
    taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
    tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
    sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA),
    Age = paste('AgeBlk', Age))

# Values of movement estimates
move_plot = ggplot(move_all_df, aes(x = To, y = value, color = model, group = model)) +
  geom_line(lwd = 1, alpha = 0.85) +
  scale_color_manual(values = colors) +
  facet_grid(Age~From) +
  theme_bw(base_size = 15) +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'To', y = "P(Movement)", color = 'Model') 

# CV of movemente estimates
move_cv_plot = ggplot(move_all_df, aes(x = To, y = cv, color = model, group = model)) +
  geom_line(lwd = 1, alpha = 0.85) +
  scale_color_manual(values = colors) +
  facet_grid(Age~From) +
  theme_bw(base_size = 15) +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'To', y = "CV P(Movement)", color = 'Model') 

ggsave(
  move_plot,
  filename = here("figs", "Manuscript_Plots", "Movement_5Area_Comparison.png"),
  width = 10, height = 6
) 

ggsave(
  move_cv_plot,
  filename = here("figs", "Manuscript_Plots", "Movement_CV_5Area_Comparison.png"),
  width = 10, height = 6
) 

### Model Diagnostics -------------------------------------------------------
##### Tag Recovery ------------------------------------------------------------
obs_tag_recovery = reshape2::melt(model_list[[1]]$obs_tag_recovery) %>% 
  rename(age = Var1, release_event = Var2, region = Var3, year = Var4, obs = value) # get observed

pred_tag_recovery = data.frame()
for(i in 1:length(model_list)) {
  if(length(dim(model_list[[i]]$obs_tag_recovery)) == 3) { # no age blocks
    pred_tag_recovery_tmp = reshape2::melt(model_list[[i]]$pred_tag_recovery) %>% 
      rename(release_event = Var1, region = Var2,  year = Var3, pred = value) %>% 
      mutate(model = model_name[i], nll = model_list[[i]]$nll[8], age = NA) # arbitrary age set for left joining
  }
  
  if(length(dim(model_list[[i]]$obs_tag_recovery)) == 4) { # includes age blocks
    pred_tag_recovery_tmp = reshape2::melt(model_list[[i]]$pred_tag_recovery) %>% 
      rename(age = Var1, release_event = Var2, region = Var3,  year = Var4, pred = value) %>% 
      mutate(model = model_name[i], nll = model_list[[i]]$nll[8])
  }
  pred_tag_recovery = rbind(pred_tag_recovery, pred_tag_recovery_tmp)
} # end i

# Join all datasets together
pred_tag_recovery = pred_tag_recovery %>% 
  group_by(year, model, region) %>% 
  summarize(pred = sum(pred, na.rm = TRUE), nll = mean(nll)) # summarize predicted

obs_tag_recovery = obs_tag_recovery %>% 
  group_by(year, region) %>% 
  summarize(obs = sum(obs, na.rm = TRUE)) # summarize observed

# join together
tag_recovery_df = pred_tag_recovery %>% 
  left_join(obs_tag_recovery, by = c("region", "year")) %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    model = factor(model, levels = model_name)
  )

# Quick relevling
tag_recovery_df = tag_recovery_df %>% 
  mutate(model = factor(model, levels = model_levels),
         # Creating facets to combine comparisons
         movement = ifelse(model %in% c("01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
         taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
         tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
         sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA)
  ) %>%  filter(model != '00-NoTag')


# Tag Recovery from 5-Area Model (movement)
movement_tag_comp = ggplot(tag_recovery_df %>% filter(!is.na(movement))) +
  geom_point(aes(x = year + 1978, y = obs), size = 2) +
  geom_line(aes(x = year + 1978, y = pred, color = model, lty= model), lwd = 1, alpha = 1) +
  coord_cartesian(ylim = c(0,500)) +
  facet_grid(movement~region) +
  scale_color_manual(values = colors[-1]) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = '', fill = '', lty = '') 

# Tag Recovery from 5-Area Model (tag likelihood)
taglike_tag_comp = ggplot(tag_recovery_df %>% filter(!is.na(taglike))) +
  geom_point(aes(x = year + 1978, y = obs), size = 2) +
  geom_line(aes(x = year + 1978, y = pred, color = model, lty= model), lwd = 1, alpha = 1) +
  coord_cartesian(ylim = c(0,500)) +
  facet_grid(taglike~region) +
  scale_color_manual(values = colors[-c(1,2,3)]) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = '', fill = '', lty = '') 

# Tag Recovery from 5-Area Model (tagrep)
tagrep_tag_comp = ggplot(tag_recovery_df %>% filter(!is.na(tagrep))) +
  geom_point(aes(x = year + 1978, y = obs), size = 2) +
  geom_line(aes(x = year + 1978, y = pred, color = model, lty= model), lwd = 1, alpha = 1) +
  coord_cartesian(ylim = c(0,500)) +
  facet_grid(tagrep~region) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = '', fill = '', lty = '') 

# Tag Recovery from 5-Area Model (sptq)
sptq_tag_comp = ggplot(tag_recovery_df %>% filter(!is.na(sptq))) +
  geom_point(aes(x = year + 1978, y = obs), size = 2) +
  geom_line(aes(x = year + 1978, y = pred, color = model, lty= model), lwd = 1, alpha = 1) +
  coord_cartesian(ylim = c(0,500)) +
  facet_grid(sptq~region) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.055, 0.625),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  labs(x = 'Year', y = 'Tag Recoveries', color = '', fill = '', lty = '') 

ggsave(
  plot_grid(movement_tag_comp, taglike_tag_comp, 
            tagrep_tag_comp, sptq_tag_comp, ncol = 1,   
            labels = c('A', 'B', 'C', 'D'),
            label_size = 18),
  filename = here("figs", "Manuscript_Plots", "TagRec_5Area_Comparison.png"),
  width = 15, height = 10
) 


##### Index -------------------------------------------------------------------
# Read in abundance index used
design_survey_index = readRDS(file = here("Data", "Survey", "regional_abundance_estimates.RDS"))
design_survey_index$area_lab = design_survey_index$NPFMC.Sablefish.Management.Area
design_survey_index = design_survey_index %>% mutate(area_lab = 
                                                       case_when(area_lab == "Aleutians" ~ "AI",
                                                                 area_lab == "Bering Sea" ~ "BS",
                                                                 area_lab == "Western Gulf of Alaska" ~ "WGOA",
                                                                 area_lab == "Central Gulf of Alaska" ~ "CGOA",
                                                                 area_lab == "Eastern Gulf of Alaska" ~ "EGOA",
                                                                 TRUE ~ area_lab)) %>% filter(Year != 2022)

# munge design based survey index
design_survey_index = design_survey_index %>% group_by(area_lab, Country, Year) %>% 
  summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), 
            se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, 
            UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)

# Get predicted idx
pred_srv = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_tmp = reshape2::melt(model_list[[i]]$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, pred = value) %>%  
    mutate(model = model_name[i], nll = model_list[[i]]$nll[5])
  pred_srv = rbind(pred_srv, pred_srv_tmp)
} # end i

# residual munging
pred_srv = pred_srv %>% 
  mutate(Year = Year + 1959, 
         area_lab = case_when(
           area_lab == 1 ~ "BS",
           area_lab == 2 ~ "AI",
           area_lab == 3 ~ "WGOA",
           area_lab == 4 ~ "CGOA",
           area_lab == 5 ~ "EGOA"),
         Country = ifelse(Country == 1, "Japan", "United States")) %>% 
  filter(Year >= 1979) 

# Quick relevling
pred_srv = pred_srv %>% 
  mutate( # Creating facets to combine comparisons
         movement = ifelse(model %in% c("00-NoTag", "01-Base", "01-Age2Move", "01-Age3Move"), 'Movement', NA),
         taglike = ifelse(model %in% c("01-Age3Move", "02-NegBin"), 'Tag Likelihood', NA),
         tagrep = ifelse(model %in% c("02-NegBin", "03-Decadal", "03-FishBlock", "03-Space"), 'Tag Reporting', NA),
         sptq = ifelse(model %in% c("03-FishBlock", "04-SptQ"), 'Spatial Q', NA)
  )


# Join datasets together
srv_bio_df = design_survey_index %>%
  left_join(pred_srv, by = c("Country", "area_lab", "Year")) %>% 
  mutate(resids = (log(sum_estimates) - log(pred)) / se,
         model = factor(model, levels = model_levels),
         area_lab = factor(area_lab, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

# Fits to index directly (movement)
movement_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(movement),Country == 'United States')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors) +
  facet_grid(movement~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.07, 0.625),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (taglike)
taglike_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(taglike),Country == 'United States')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1, 2,3)]) +
  facet_grid(taglike~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.07, 0.625),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (tagrep)
tagrep_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(tagrep),Country == 'United States')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.07, 0.625),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (sptq)
sptq_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(sptq),Country == 'United States')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.07, 0.625),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  plot_grid(movement_idx_plot, taglike_idx_plot, tagrep_idx_plot, sptq_idx_plot, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 23),
  filename = here("figs", "Manuscript_Plots", "IdxFits_US_5Area_Comparison.png"),
  width = 23, height = 20
) 

# Japan fits
# Fits to index directly (movement)
movement_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(movement),Country == 'Japan')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors) +
  facet_grid(movement~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.93, 0.85),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (taglike)
taglike_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(taglike),Country == 'Japan')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1, 2,3)]) +
  facet_grid(taglike~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.93, 0.85),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (tagrep)
tagrep_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(tagrep),Country == 'Japan')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1,2,3,4)]) +
  facet_grid(tagrep~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.93, 0.85),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

# Fits to index directly (sptq)
sptq_idx_plot = ggplot(srv_bio_df %>% filter(!is.na(sptq),Country == 'Japan')) +
  geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI), alpha = 0.6, size = 1) +
  geom_line(aes(x = Year, y = pred, color = model, lty = model), alpha = 1, lwd = 2) +
  scale_color_manual(values = colors[-c(1,2,3,4,5,6,8)]) +
  facet_grid(sptq~area_lab, scales = "free") +
  labs(x = "Year", y = "Relative Population Numbers", color = "", lty = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.93, 0.85),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, -0.075, 0.1), "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  plot_grid(movement_idx_plot, taglike_idx_plot, tagrep_idx_plot, sptq_idx_plot, ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 23),
  filename = here("figs", "Manuscript_Plots", "IdxFits_JP_5Area_Comparison.png"),
  width = 23, height = 20
) 

