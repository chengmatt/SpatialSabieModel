# Purpose: To plot 5 Area 1960 models for comparisons when building complexity in the model
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date: 9/6/24

# Set up ------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
library(here)
source(here("R", "Utility_Fxns.R"))
source(here("R", "Model_Compare_Fxns.R"))

out_path = here("Output", "Final Models", "5-Area-1960")

# Look at convergence real quick
files = list.files(out_path)
files = files[!str_detect(files, ".RDS")]
conv_all = data.frame()
for(i in 1:length(files)) {
  mod = readRDS(here(out_path, files[i], 'sd_report.RDS'))
  tmp = data.frame(mod = files[i], pdHess = mod$pdHess, grad = max(abs(mod$gradient.fixed)))
  conv_all = rbind(tmp ,conv_all)
}

# 01 - Movement (Age vs. Time) -----------------------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-Base"),
    here(out_path, "5-Area-1960-01-Age_2_Move"),
    here(out_path, "5-Area-1960-01-Age_3_Move"),
    here(out_path, "5-Area-1960-01-Time_2_Move"),
    here(out_path, "5-Area-1960-01-Time_4_Move")
  ),
  model_name = c("01-Base", "01-Age2Move", "01-Age3Move", "01-Time2Move", "01-Time4Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01-Age_vs_TimeMove.pdf"
)


# 01 - Movement (Age x Time) ----------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-Age_3_Move"),
    here(out_path, "5-Area-1960-01-Age_2_Time_2_Move"),
    here(out_path, "5-Area-1960-01-Age_3_Time_2_Move"),
    here(out_path, "5-Area-1960-01-Age_2_Time_3_Move")
  ),
  model_name = c("01-Age3Move", "01-Age2Time2Move", "01-Age3Time2Move", "01-Age_2_Time_3_Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01-Age_vs_AgexTime.pdf"
)

# 02 - Likelihoods ----------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-Age_3_Move"),
    here(out_path, "5-Area-1960-02-NegBin"),
    here(out_path, "5-Area-1960-02-Mltnml")
  ),
  model_name = c("01-Age3Move (Pois)", "02-NegBin", "02-Mltnml"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "02-Likelihoods.pdf"
)


# 03 - Reporting Rates ----------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-02-NegBin"),
    here(out_path, "5-Area-1960-03-Decadal"),
    here(out_path, "5-Area-1960-03-FishBlock"),
    here(out_path, "5-Area-1960-03-Space")
  ),
  model_name = c("02-NegBin (Const)", "03-Decadal", "03-FishBlock", "03-Space"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "03-RepRates.pdf"
)

# 04 - Spatial Q ----------------------------------------------

# Spatial Q changes the scaling of stuff. In a spatially invariant q case,
# scaling is based on realtive difference in CPUE. In a spatially varying q case,
# scaling is instead based on a combination of the index and catch information primarily.
# So the model will attempt to estimate a q that fits the index well and then maybe
# 'overleverage' the catch information. By contrast, a spatially invariant q
# will probably fit index less well, but overlerverage the relative CPUE scale probably.

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-03-FishBlock"),
    here(out_path, "5-Area-1960-04-SptQ")
  ),
  model_name = c("03-FishBlock", "04-SptQ"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "04-SptQ.pdf"
)


# 00 - No Tag Scenario ----------------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-03-FishBlock"),
    here(out_path, "5-Area-1960-00-NoTag")
  ),
  model_name = c("03-FishBlock", "00-NoTag"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "00-NoTag.pdf"
)


# Compare Final to 1area --------------------------------------------------
area1_dat = readRDS(here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final", "data.RDS"))
area1 = readRDS(here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final", "mle_report.RDS"))
area5_dat = readRDS(here(out_path, "5-Area-1960-03-FishBlock", "data.RDS"))
area5 = readRDS(here(out_path, "5-Area-1960-04-SptQ", "mle_report.RDS"))

par(mfrow = c(1,3))
# SSB
plot(1960:2021, area1$SSB_yr, ylim = c(0, 300), type = 'l', lwd = 3, xlab = "Year", ylab = "SSB")
lines(1960:2021, rowSums(area5$SSB_yr), ylim = c(0, 300), type = 'l', lwd = 3, col = "blue")
legend(x = 1990, y = 90, legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Depletion
plot(1960:2021, area1$SSB_yr / area1$SSB_yr[1], type = 'l', lwd = 3, xlab = "Year", ylab = "SSB", ylim = c(0, 1))
lines(1960:2021, rowSums(area5$SSB_yr) / rowSums(area5$SSB_yr)[1],  type = 'l', lwd = 3, col = "blue")
legend(x = 1990, y = 0.2, legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Rec
plot(1960:2021, area1$recruitment_yr, ylim = c(0, 200), type = 'l', lwd = 3, xlab = "Year", ylab = "Rec")
lines(1960:2021, rowSums(area5$recruitment_yr), ylim = c(0, 200), type = 'l', lwd = 3, col = "blue")
legend(x = 1990, y = 90, legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Harvest rate

# Get weight at age to calcualte harvest rate
waa_f = array(area1_dat$female_mean_weight_by_age, dim = c(length(area1_dat$ages), area1_dat$n_regions, length(area1_dat$years) + 1))
waa_m = array(area1_dat$male_mean_weight_by_age, dim = c(length(area1_dat$ages), area1_dat$n_regions, length(area1_dat$years) + 1))

# Get array
# females
sel_fixed_f = aperm(replicate(1, cbind(replicate(length(1960:2015), area1$sel_fixed_f[,1]),
                                       replicate(length(2016:2022), area1$sel_fixed_f[,2]))), c(1,3,2))

# males
sel_fixed_m = aperm(replicate(1, cbind(replicate(length(1960:2015), area1$sel_fixed_m[,1]),
                                       replicate(length(2016:2022), area1$sel_fixed_m[,2]))), c(1,3,2))

# 1 Area harvest rate
# Get exploitable biomass for fixed gear at age
expl_fixed_age1 = sel_fixed_f * area1$natage_f * waa_f +
                  sel_fixed_m * area1$natage_m * waa_m 
pred_fixed_cat1 = area1$fixed_fishery_catch # get predicted fixed gear catch
# Get exploitable biomass for trawl gear at age
expl_trwl_age1 = as.vector(area1$sel_trwl_f) * area1$natage_f * waa_f +
  as.vector(area1$sel_trwl_m) * area1$natage_m * waa_m 
pred_trawl_cat1 = area1$trwl_fishery_catch # get predicted fixed gear catch
# Sum across ages
expl_fixed1 = apply(expl_fixed_age1, c(2,3), sum)
expl_trwl1 = apply(expl_trwl_age1, c(2,3), sum)

# 5 Area harvest rate
waa_f = array(area5_dat$female_mean_weight_by_age, dim = c(length(area1_dat$ages), area5_dat$n_regions, length(area5_dat$years) + 1))
waa_m = array(area5_dat$male_mean_weight_by_age, dim = c(length(area1_dat$ages), area5_dat$n_regions, length(area5_dat$years) + 1))

# females
sel_fixed_f = aperm(replicate(5, cbind(replicate(length(1960:2015), area5$sel_fixed_f[,1]),
                                       replicate(length(2016:2022), area5$sel_fixed_f[,2]))), c(1,3,2))

# males
sel_fixed_m = aperm(replicate(5, cbind(replicate(length(1960:2015), area5$sel_fixed_m[,1]),
                                       replicate(length(2016:2022), area5$sel_fixed_m[,2]))), c(1,3,2))

expl_fixed_age5 = sel_fixed_f * area5$natage_f * waa_f +
                  sel_fixed_m * area5$natage_m * waa_m 
pred_fixed_cat5 = area5$fixed_fishery_catch # get predicted fixed gear catch
# Get exploitable biomass for trawl gear at age
expl_trwl_age5 = as.vector(area5$sel_trwl_f) * area5$natage_f * waa_f +
  as.vector(area5$sel_trwl_m) * area5$natage_m * waa_m 
pred_trawl_cat5 = area5$trwl_fishery_catch # get predicted fixed gear catch
# Sum across ages
expl_fixed5 = apply(expl_fixed_age5, c(2,3), sum)
expl_trwl5 = apply(expl_trwl_age5, c(2,3), sum)

plot(1960:2021, t(pred_fixed_cat1 / expl_fixed1[-63]), col = 'black', type = 'l')
lines(1960:2021, colSums(pred_fixed_cat5) / colSums(expl_fixed5[,-63]), col = 'blue')
legend(x = 1990, y = 0.08, legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

