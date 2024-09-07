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

# 01-Poisson vs. 01-NegBin ------------------------------------------------

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-Poisson"),
    model_2 = here(out_path, "5-Area-1960-01-NegBin")
  ),
  model_name = c("01-Pois", "01-NegBin"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_Pois_vs_NegBin.pdf"
)

# Negative Binomial seems to fit index data better, at the cost of fitting recapture data (makes sense)

# 01-NegBin (Constant Tag Recovery) vs. 02-Recent Tag Recovery ------------------------------------------------

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-02-Recent")
  ),
  model_name = c("01-NegBin-Constant", "02-Recent"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_02-Recent.pdf"
)

# Comparing constant to recent time-blocked tag recovery, there does not seem
# to be a huge difference in either fits to indices, tag recapture data, and
# general dynamics seem similar. Doesn't seme necessary to incorporate

# Circle back to this (where is information from tag recovery coming from??)

# 01-NegBin (Constant Tag Recovery) vs. 02-Decadal Tag Recovery ------------------------------------------------

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-02-Decadal")
  ),
  model_name = c("01-NegBin-Constant", "02-Decadal"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_02-Decadal.pdf"
)

# Again, not much difference with time-varying vs. time-constant tag recovery rates

# 01-NegBin (Constant Tag Recovery) vs. 02-Spatial Tag Recovery ------------------------------------------------
compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-02-Space")
  ),
  model_name = c("01-NegBin-Constant", "02-Space"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_02-Space.pdf"
)

# Again, not much difference with time-varying vs. constant tag recovery rates
# Seems like we can continue assuming constant spatial and time rates, without large consequences

# 03-NegBin (Constant Tag Recovery) vs. 03-Time Block (4) Movement  ------------------------------------------------
movement_blk_4_rep <- readRDS(here(out_path, "5-Area-1960-03-Time_4_Move", "mle_report.RDS"))
movement_blk_4_sd <- readRDS(here(out_path, "5-Area-1960-03-Time_4_Move", "sd_report.RDS"))

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-03-Time_4_Move")
  ),
  model_name = c("01-NegBin-Constant", "03-Time4Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_03_Time4Move.pdf"
)

# Movement parameters for 4 time blocks is pretty uncertain, but fits
# to indices and tag recapture data improve slightly. Given high uncertainty
# in movement dynamics, it may not be great to proceed with 4 time blocks (look at alternative blocking for time)

# 01-NegBin (Constant Tag Recovery) vs. 03-Time Block (3) Movement  ------------------------------------------------
movement_blk_3_rep <- readRDS(here(out_path, "5-Area-1960-03-Time_3_Move", "mle_report.RDS"))
movement_blk_3_sd <- readRDS(here(out_path, "5-Area-1960-03-Time_3_Move", "sd_report.RDS"))

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-03-Time_3_Move")
  ),
  model_name = c("01-NegBin-Constant", "03-Time3Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_03_Time3Move.pdf"
)

# Fits to tag recapture data improve, as well as to indices. Movement parameters are
# less uncertain, and this is a more flexible parmaeterization than 2 time blocks,
# Likely will go forward with this in future investigations.

# 01-NegBin (Constant Tag Recovery) vs. 03-Time Block (2) Movement  ------------------------------------------------
movement_blk_2_rep <- readRDS(here(out_path, "5-Area-1960-03-Time_2_Move", "mle_report.RDS"))
movement_blk_2_sd <- readRDS(here(out_path, "5-Area-1960-03-Time_2_Move", "sd_report.RDS"))

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-03-Time_2_Move")
  ),
  model_name = c("01-NegBin-Constant", "03-Time2Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBin_vs_03_Time2Move.pdf"
)

# Fits to tag recapture data improve ever so slightly but not much improvement, 
# neither for indices. Movement dynmiacs are less uncertain and they seem fairly reliable.

