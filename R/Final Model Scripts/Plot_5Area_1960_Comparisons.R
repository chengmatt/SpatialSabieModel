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

# Decadal tag recovery fits to tag recaptures better (given increased flexibility), but fits to other
# data sources, as well as key quantities did not change substatially. Does not seem worth the extra complexity

# 03-Movement 3 Time Blocks vs. 03-Movement 4 Time Blocks ------------------------------------------------

# Looking at these two models, the movement parameters are so uncertain that 
# I would argue they're unreliable and likely on a bound, and cannot be 
# estimated adequately (hence, not even going to 
# look at the plots). However, the 2 time block models seem to be estimated relatively well. 

# 01-NegBin (Constant) vs. 03-Movement 2 Time Blocks ------------------------------------------------

compare_2_models(
  model_path = list(
    model_1 = here(out_path, "5-Area-1960-01-NegBin"),
    model_2 = here(out_path, "5-Area-1960-03-Time_2_Move")
  ),
  model_name = c("01-NegBin (ConstMove)", "03-Time2Move"),
  fig_path = here("Figs", "5-Area-1960"),
  fig_name = "01_NegBinConst_vs_Time2Move.pdf"
)

# Comparing a time-blocked movement model against a time-invariant movement model,
# it seems like fits to indices improve a decent bit, as well as fits to recapture data.
# However, there seems to be major differences in regional recruitment and scaling SSB (although
# the absolute scale of SSB seems similar)

# Need to come back to this and figure out why its doing that...
 

model_1 = readRDS(here(out_path, "5-Area-1960-01-NegBin", "mle_report.RDS"))
model_2 = readRDS(here(out_path, "5-Area-1960-03-Time_2_Move", "mle_report.RDS"))

reshape2::melt(model_1$movement_matrix) %>%
  rename(From = Var1, To = Var2, TimeBlock = Var3) %>%
  mutate(
    From = case_when(
      From == 1 ~ "BS",
      From == 2 ~ "AI",
      From == 3 ~ "WGOA",
      From == 4 ~ "CGOA",
      From == 5 ~ "EGOA"
    ), From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), To = factor(To, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  ggplot(aes(x = From, y = To, fill = value, label = round(value, 3))) +
  geom_tile(alpha = 0.8) +
  geom_text(size = 5) +
  scale_fill_viridis_c() +
  facet_wrap(~TimeBlock) +
  labs(title = "Constant Move")

reshape2::melt(model_2$movement_matrix) %>%
  rename(From = Var1, To = Var2, TimeBlock = Var3) %>%
  mutate(
    From = case_when(
      From == 1 ~ "BS",
      From == 2 ~ "AI",
      From == 3 ~ "WGOA",
      From == 4 ~ "CGOA",
      From == 5 ~ "EGOA"
    ), From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), To = factor(To, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  ggplot(aes(x = From, y = To, fill = value, label = round(value, 3))) +
  geom_tile(alpha = 0.8) +
  geom_text(size = 5) +
  scale_fill_viridis_c() +
  facet_wrap(~TimeBlock) +
  labs(title = "Time Block Move")
