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
  mod = readRDS(here(out, files[i], 'sd_report.RDS'))
  tmp = data.frame(mod = files[i], pd = mod$pdHess, grad = max(abs(mod$gradient.fixed)))
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
    here(out_path, "5-Area-1960-01-Age_3_Time_2_Move")
  ),
  model_name = c("01-Age3Move", "01-Age2Time2Move", "01-Age3Time2Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01-Age_vs_AgexTime.pdf"
)

# 03 - Reporting Rates ----------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-01-Age_3_Move"),
    here(out_path, "5-Area-1960-03-Decadal"),
    here(out_path, "5-Area-1960-03-FishBlock"),
    here(out_path, "5-Area-1960-03-Space")
  ),
  model_name = c("01-Age3Move (Const)", "03-Decadal", "03-FishBlock", "03-Space"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "03-RepRates.pdf"
)

# 04 - Spatial Q ----------------------------------------------

compare_5sptl_models(
  model_path = list(
    here(out_path, "5-Area-1960-03-FishBlock"),
    here(out_path, "5-Area-1960-04-SptQ")
  ),
  model_name = c("03-FishBlock", "04-SptQ"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "04-SptQ.pdf"
)

