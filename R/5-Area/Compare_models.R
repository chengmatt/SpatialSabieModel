#'
#' Compare all 5 area models using bookdown code
#'
#'
source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)

###########
## Repeat for the direct-ageing AF inputs
##########
model_description = '
- "Init"  The initial model, assumes multinomial likelihoods for composition. AFs were estimated using the Age length key. Tag recovery data was assumed to be Negative Binomial. Annual recruitment deviations are shared across the regions
- "D-M" The same "Init" but assumes compositional data is Dirichlet Multinomial
- "Regional Recruitment" same as "D-M" but with regional specific annual recruitment deviations
- "Constrained recruitment" same as "Regional Recruitment" but recruit devs in each region are constrained to sum = 0
- "tag-poisson" same as "Constrained recruitment" but tag-likelihood assumed to be Poisson instead of Negative Binomial
- "decadal-RR" same as "tag-poisson" but tag-reporting rate 
- "decadal-RR-NB" same as "decadal-RR" but changed back to Negative binomial tag-likelihood
- "CHG trwl sel" Change the trawl selectivity to be a 3 parameter double normal

'
bookdown_labels = c("Init","D-M",  "Regional recruitment", "Constrained recruitment", "tag-poisson", "decadal-RR", "decadal-RR-NB", "CHG trwl sel")

model_labels = c("Model_01_srva", "Model_02_srva", "Model_03_srva", "Model_04_srva", "Model_05_srva", "Model_06_srva", "Model_07_srva", "Model_08_srva")
rbind(bookdown_labels,model_labels)

model_dir = file.path(DIR$app_5A,model_labels)

bookdown_dir = file.path(DIR$book5A, "Comparison_individual_DirectAgeing")

summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                            model_description = model_description)
bookdown_dir = file.path(DIR$book5A, "Comparison_together_DirectAgeing")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                          model_description = model_description)

## 1960 model
model_labels = c("Model1960_01_srva", "Model1960_02_srva", "Model1960_03_srva")#, "Model1960_04_srva", "Model1960_05_srva", "Model1960_06_srva")
bookdown_labels = c("1960","decade tag",  "spatial q")#, "S2", "S2 10 YCS", "S2 const_rec")

model_dir = file.path(DIR$app_5A,model_labels)
bookdown_dir = file.path(DIR$book5A, "Comparison_individual_1960")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                            model_description = model_description)
bookdown_dir = file.path(DIR$book5A, "Comparison_together_1960")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                          model_description = model_description)



bookdown_labels = c("1960", "1977", "1990")
model_labels = c( "Model1960_06a","Model_07a", "Model1990_01a")
model_dir = file.path(DIR$app_5A,model_labels)

bookdown_dir = file.path(DIR$book5A, "Comparison_together_DirectAgeing_multi_year")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                          model_description = model_description)

bookdown_dir = file.path(DIR$book5A, "Comparison_indiviudal_DirectAgeing_multi_year")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, 
                            model_description = model_description)

# options(warn=0) ## if warnings are being set to errors
fig_path = file.path(DIR$app_5A, "Comparison")
if(!dir.exists(fig_path)) {
  dir.create(fig_path)
}
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))


scenarios = c("Model_01", "Model_01a", "Model_02","Model_03", "Model_04")
labels = c("Init (ALK)","Direct\nAgeing", "D-M", "Regional\nRecruitment", "Recuit\nconstrained")
mle_ls = list()
for(i in 1:length(labels)) 
  mle_ls[[labels[i]]] = readRDS(file.path(DIR$app_5A, scenarios[i], "mle_report.RDS"))

## get derived quantities
ssbs = get_multiple_ssbs(mle_ls, labels, region_key = region_key, depletion = F)
ssbs_percent = get_multiple_ssbs(mle_ls, labels, region_key = region_key, depletion = T)
Fs = get_multiple_Fs(mle_ls, labels, region_key = region_key)
rec_df = get_multiple_recruits(mle_ls, labels, region_key = region_key)
sel_df = get_multiple_selectivities(mle_ls, labels)

## plot SSBs
ggplot(ssbs %>% group_by(Year, label), aes(x = Year, y = SSB, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  labs(x = "Year", y = "SSB (kt)", col = "Model", linetype = "Model") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
        ) +
  facet_wrap(~Region)

ggsave(filename = file.path(fig_path, "SSB_summarise_over_models.png"), width = 9, height = 7)


## plot percent SSBs
ggplot(ssbs_percent %>% group_by(Year, label), aes(x = Year, y = Depletion, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  labs(x = "Year", y = "Depletion (SSB/B0)", col = "Model", linetype = "Model") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
  ) +
  facet_wrap(~Region)

ggsave(filename = file.path(fig_path, "Depletion_summarise_over_models.png"), width = 9, height = 7)

## plot Fs
# fixed gear Fs
ggplot(Fs %>% group_by(Year, label) %>% filter(Fishery == "Fixed gear"), aes(x = Year, y = F, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  ggtitle("Fixed gear") +
  labs(x = "Year", y = "F", col = "Model", linetype = "Model") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
  ) +
  facet_wrap(~Region)

ggsave(filename = file.path(fig_path, "Fixed_gear_Fs.png"), width = 9, height = 7)

# Trawl gear Fs
ggplot(Fs %>% group_by(Year, label) %>% filter(Fishery == "Trawl"), aes(x = Year, y = F, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  ggtitle("Trawl gear") +
  labs(x = "Year", y = "F", col = "Model", linetype = "Model") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
  ) +
  facet_wrap(~Region)

ggsave(filename = file.path(fig_path, "Trawl_gear_Fs.png"), width = 9, height = 7)

## recruits
ggplot(rec_df %>% group_by(Year, label), aes(x = Year, y = Recruitment, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ylim(0,NA) +
  labs(x = "Year", y = "Recruits (million)", col = "Model", linetype = "Model") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
  ) +
  facet_wrap(~Region)

ggsave(filename = file.path(fig_path, "Recruitment.png"), width = 9, height = 7)


## Selectivities
ggplot(sel_df, aes(x = age, y = value, col = label, linetype = label)) +
       geom_line(linewidth = 1.1) +
         theme_bw() +
         ylim(0,NA) +
         labs(x = "Year", y = "Selectivity", col = "Model", linetype = "Model") +
         theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 14),
               legend.position = "bottom"
         ) +
         facet_wrap(~name)
ggsave(filename = file.path(fig_path, "Selectivities.png"), width = 9, height = 7)

       