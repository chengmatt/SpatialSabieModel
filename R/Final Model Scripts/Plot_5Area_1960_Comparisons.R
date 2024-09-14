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

# 01-Tag Likelihoods ------------------------------------------------

compare_sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-Poisson"),
    here(out_path, "5-Area-1960-01-NegBin"),
    here(out_path, "5-Area-1960-01-Multinomial")
  ),
  model_name = c("01-Pois", "01-NegBin", "01-Mltnml"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_TagLikelihoods.pdf"
)


# 02 - Reporting Rates ----------------------------------------------------

compare_sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-NegBin"),
    here(out_path, "5-Area-1960-02-Recent"),
    here(out_path, "5-Area-1960-02-FishBlock"),
    here(out_path, "5-Area-1960-02-Space"),
    here(out_path, "5-Area-1960-02-Decadal")
  ),
  model_name = c("01-NegBin (Constant)", "02-Recent", "02-FishBlock", "02-Spatial", "02-Decadal"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "02_TagReporting.pdf"
)

# 03 - Movement ----------------------------------------------------

compare_sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-02-FishBlock"),
    here(out_path, "5-Area-1960-03-Time_4_Move"),
    here(out_path, "5-Area-1960-03-Time_3_Move"),
    here(out_path, "5-Area-1960-03-Time_2_Move")
  ),
  model_name = c("02-FishBlock (Constant)", "03-Move4", "03-Move3", "03-Move2"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "03_Movement.pdf"
)


# 04 - Spatial Q ----------------------------------------------------

compare_sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-03-Time_2_Move"),
    here(out_path, "5-Area-1960-04-SptQ")
  ),
  model_name = c("03-Move2", "04-SptQ"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "04_SptQ.pdf"
)


# Time Series - 1Area vs 5Area --------------------------------------------

area_1 = readRDS(here("Output", "Final Models", "1-Area-1960", "1-Area-1960-Final", "mle_report.RDS"))
area_5 = readRDS(here(out_path, "5-Area-1960-04-SptQ", "mle_report.RDS"))

par(mfrow = c(1,2))
# SSB
plot(1960:2021, area_1$SSB_yr, ylim = c(0, 315), type = 'l', xlab = "Year", ylab = "SSB", lwd = 3)
lines(1960:2021, rowSums(area_5$SSB_yr), type = 'l', xlab = "Year", ylab = "SSB", col = "blue", lwd = 3)
legend(x = c(1960, 1960), y = c(50, 50), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Recruitment
plot(1960:2021, area_1$recruitment_yr, ylim = c(0, 150), type = 'l', xlab = "Year", ylab = "Recruitment", lwd = 3)
lines(1960:2021, rowSums(area_5$recruitment_yr), type = 'l', xlab = "Year", ylab = "SSB", col = "blue", lwd = 3)
legend(x = c(1990, 1990), y = c(120, 120), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

par(mfrow = c(2,3))
# F Fixed Gear
plot(1960:2021, colSums(area_5$annual_F_fixed), type = 'l',col = "blue", lwd = 3,  xlab = "Year", ylab = "F (Fixed Gear)")
lines(1960:2021, area_1$annual_F_fixed,  type = 'l', lwd = 3)
legend(x = c(1980, 1980), y = c(0.8, 0.8), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Proxy Fixed Gear Harvest Rate
plot(1960:2021, t(area_1$fixed_fishery_catch) / area_1$total_biomass_yr, type = 'l', col = "black", lwd = 3, xlab = "Year", ylab = "Proxy Harvest Rate (Fixed Gear - Catch / Total Biomass)")
lines(1960:2021, rowSums(t(area_5$fixed_fishery_catch)) / rowSums(area_5$total_biomass_yr), col = "blue", lwd = 3)
legend(x = c(1980, 1980), y = c(0.02, 0.02), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Correlation between F spatial vs. non spatial
plot(colSums(area_5$annual_F_fixed), area_1$annual_F_fixed, xlab = "F (Fixed Gear Spatial)", ylab = "F (Fixed Gear 1Area)")

# F Trawl Gear
plot(1960:2021, colSums(area_5$annual_F_trwl), type = 'l',col = "blue", lwd = 3,  xlab = "Year", ylab = "F (Trawl Gear)")
lines(1960:2021, area_1$annual_F_trwl,  type = 'l', lwd = 3)
legend(x = c(1980, 1980), y = c(0.8, 0.8), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Proxy Fixed Gear Harvest Rate
plot(1960:2021, t(area_1$annual_trwl_catch_pred) / area_1$total_biomass_yr, type = 'l', col = "black", lwd = 3, xlab = "Year", ylab = "Proxy Harvest Rate (Trawl Gear - Catch / Total Biomass)")
lines(1960:2021, rowSums(t(area_5$annual_trwl_catch_pred)) / rowSums(area_5$total_biomass_yr), col = "blue", lwd = 3)
legend(x = c(1980, 1980), y = c(0.04, 0.04), legend = c("1-Area", "5-Area"), fill = c("black", "blue"))

# Correlation between F spatial vs. non spatial
plot(colSums(area_5$annual_F_trwl), area_1$annual_F_trwl, xlab = "F (Trawl Gear Spatial)", ylab = "F (Trawl Gear 1Area)")
