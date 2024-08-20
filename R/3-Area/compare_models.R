#'
#' Compare all 3 area models using bookdown code
#'
#'
source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
## 1977 models
model_labels = c("NoTagData_countrysrv_regionaldevsa", "NoTagData_est_move_rec_srva", "fenske_movement_countrysrva",  "Model_05_srva", "Model_06_srva")#, "Model_02a")#, "1960_TagData_01a", "1960_TagData_02a","1960_TagData_03a")
bookdown_labels = c("reg rec", "est move", "fenske move", "reduc rec", "spatial q")#, "NB reg rec")#, "Poisson tag","NB tag", "decade TR", "reduced_initage)
model_description = '
- "reg rec"  Same as "1977" but estimate recruitment specific deviations
- "est move"  Estimate movement with no tagging data"
- "fenske move"  Fix movement to estimates from Fenske 2023"
- "Poisson" Same as "est move" but includeds tag-data assuming a Poisson likelihood
- "Pois rec" Same as "Poisson" but regional recruitment
- "dec TR" Decadal teg-reporting
- "NB" Change tag-likelihood to negative binomial
- "reduc rec" only estimate the first 15 inital age-devs
'
model_dir = file.path(DIR$app_3A,model_labels)
bookdown_dir = file.path(DIR$book3A, "Comparison_individual_1977_init")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book3A, "Comparison_together_1977_init")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
## 1960 models
model_labels = c("1960_NoTagData_regrecruit_srva", "1960_NoTagData_est_move_srva","1960_fenske_movement_countrysrva","1960_TagData_01_srva", "1960_TagData_03_srva", "1960_spatial_qsa")#, "Model_02a")#, "1960_TagData_01a", "1960_TagData_02a","1960_TagData_03a")
bookdown_labels = c("reg rec","est_move", "fenske","Poisson","NB", "spatialq")#, "NB reg rec")#, "Poisson tag","NB tag", "decade TR", "reduced_initage)
model_description = '
- "1960"  no movement, no tagging data, comp likelihood assumed multinomial, global rec devs,
- "reg rec"  Same as "1960" but estimate recruitment specific deviations
- "est move"  Estimate movement with no tagging data"
- "fenske move"  Fix movement to estimates from Fenske 2023"
- "Poisson" Same as "est move" but includeds tag-data assuming a Poisson likelihood
- "NB" Negative Binomial
'
model_dir = file.path(DIR$app_3A,model_labels)
bookdown_dir = file.path(DIR$book3A, "Comparison_individual_1960_init")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book3A, "Comparison_together_1960_init")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

## pick the best 1977 to compare with the best 1960 model
model_labels = c("Model_05_srva", "1960_TagData_03_srva", "1960_fenske_movement_spatial_q_countrysrva")
bookdown_labels = c("1977", "1960")#, "NB reg rec")#, "Poisson tag","NB tag", "decade TR", "reduced_initage)

model_description = '
- "1977"  Three area model starting in 1977
- "1960"  Three area model starting in 1960
- "Fenske"  Three area model starting in 1960 using fixed movement from Fenske
'
model_dir = file.path(DIR$app_3A,model_labels)
bookdown_dir = file.path(DIR$book3A, "Comparison_individual_1960_vs_1977")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book3A, "Comparison_together_1960_vs_1977")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

## bring in Fenske's SSBs and compare
fenske_ssb = readxl::read_xlsx(file.path(DIR$data, "FromFenske", "ssb data for craig.xlsx"), sheet = 2)
colnames(fenske_ssb) = c("Year", "BS_AI_WGOA", "CGOA", "EGOA")
fenske_srv = read.csv(file.path(DIR$data, "FromFenske", "japanese_srv_inputs.csv"))

mod_1977 = readRDS(file.path(DIR$app_3A, model_labels[1], "mle_report.RDS"))
data = readRDS(file.path(DIR$app_3A, model_labels[1], "data.RDS"))
mod_1960 = readRDS(file.path(DIR$app_3A, model_labels[2], "mle_report.RDS"))
region_key = readRDS(file.path(DIR$app_3A, model_labels[2], "region_key.RDS"))

ssb_77 = get_SSB(mod_1977, region_key = region_key)
ssb_77$run = "1977"
ssb_60 = get_SSB(mod_1960, region_key = region_key)
ssb_60$run = "1960"
ssb_fenske = fenske_ssb %>% pivot_longer(!Year)
colnames(ssb_fenske) = c("Year", "Region", "SSB")
ssb_fenske$run = "Fenske"


plot_index_fits(mod_1960, survey_labels = c("Japanese", "Domestic"))
full_ssb = rbind(ssb_77, ssb_60, ssb_fenske)

ggplot(full_ssb, aes(x = Year, y = SSB, col = run, linetype = run)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~Region) +
  theme_bw() +
  ylim(0,NA)
ggsave(filename = file.path(DIR$app_3A, "regional_ssb_comparison.png"))

plot_movement(mod_1977, region_key)
plot_movement(mod_1977, region_key)


global_ssb = full_ssb %>% group_by(Year, run) %>% summarise(value = sum(SSB))
ggplot(data = global_ssb) +
  geom_line(aes(x = Year, y = value, col = run, linetype = run), linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA)
ggsave(filename = file.path(DIR$app_3A, "global_ssb_comparison.png"))

ndx_77 = get_index(mod_1977, region_key = region_key, survey_labels = c("Japanese", "Domestic"))
ndx_60 = get_index(mod_1960, region_key = region_key, survey_labels = c("Japanese", "Domestic"))
ndx_77$run = "1977"
ndx_60$run = "1960"

fenske_ndx = fenske_srv %>% pivot_longer(!year) %>% mutate(name = 
                                                            case_when(name == "BSAIWG" ~ "BS_AI_WGOA",
                                                                      name == "CG" ~ "CGOA",
                                                                      name == "EG" ~ "EGOA",
                                                                      TRUE ~ name))
## sum
colnames(fenske_ndx) = c("Year", "Region", "Observed")
fenske_ndx$run = "Fenske"
fenske_ndx$Observed = fenske_ndx$Observed * 1000
full_ndx = rbind(fenske_ndx, ndx_60 %>% filter(Survey == "Japanese") %>% select(colnames(fenske_ndx)),  ndx_77 %>% filter(Survey == "Japanese") %>% select(colnames(fenske_ndx)))

ggplot(data = full_ndx) +
  geom_line(aes(x = Year, y = Observed, col = run, linetype = run), linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  facet_wrap(~Region)
## 
mod_1960$mean_rec
mod_1960$Bzero
mod_1977$mean_rec
mod_1977$Bzero
## why large R0 but not B0
init_age = calculate_initial_numbers_at_age(n_regions = 3, n_ages = 30, R0 = mod_1960$mean_rec / 2, movement_matrix = mod_1960$movement_matrix[,,1], natural_mortality = data$M[,1])

plot(data$maturity[,1])
lines(mod_1960$sel_srv_f[,1,1], col = "red")
lines(mod_1960$sel_srv_m[,1,1], col = "blue")

plot(data$maturity[,1])
lines(mod_1960$sel_fixed_f[,1], col = "red")
lines(mod_1960$sel_fixed_m[,1], col = "blue")

move = get_movement(mod_1960, region_key = region_key)
levels = c("BS_AI_WGOA", "CGOA", "EGOA")

move %>% pivot_wider(id_cols = From, names_from = To, values_from = Proportion)

###########
## Lets look at tag fits
###########
# mod_1977
# mod_1960
region_key$area = c("WG", "CG", "EG")
tags_60 = get_tag_recovery_obs_fitted_values(mod_1960, region_key)
tags_77 = get_tag_recovery_obs_fitted_values(mod_1977, region_key)
tags_60$residual = tags_60$observed - tags_60$predicted
tags_60$recovery_region = factor(tags_60$recovery_region, levels = region_key$area)
tags_60$release_region = factor(tags_60$release_region, levels = region_key$area)

index_fits = get_index(mod_1960, region_key = region_key, survey_labels = c("Japanese", "Domestic")) 
ggplot(index_fits, aes(x = Year)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth= 1.1, linetype = "dashed") +
  geom_errorbar(aes(ymin=L_CI, ymax=U_CI, col = "Observed"), width=.2, position=position_dodge(.9)) +
  guides( linewidth = "none", linetype = "none") +
  labs(y = "Index", col = "", linetype = "") +
  facet_grid(Survey ~ Region) +
  theme_bw() 
ggsave(filename = file.path(DIR$app_3A, "survey_ndx.png"), width = 12, height = 10)

ggplot(tags_60 %>% filter(release_year %in% 1978:1985), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(DIR$app_3A, "tag_fits_60_1.png"), width = 10, height = 10)

ggplot(tags_60 %>% filter(release_year %in% 1986:1995), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(DIR$app_3A, "tag_fits_60_2.png"), width = 10, height = 10)

ggplot(tags_60 %>% filter(release_year %in% 1996:2005), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(DIR$app_3A, "tag_fits_60_3.png"), width = 10, height = 10)

ggplot(tags_60 %>% filter(release_year %in% 2006:2015), aes(x = recovery_region, y = release_region, fill = residual)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "#FFFFCC",
                       low = "#FF0000") +
  facet_grid(release_year ~ recovery_year)
ggsave(filename = file.path(DIR$app_3A, "tag_fits_60_4.png"), width = 10, height = 10)
