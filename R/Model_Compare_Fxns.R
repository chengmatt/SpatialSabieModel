#' Title Compare 5-area spatial models
#'
#' @param model_path Path to models
#' @param model_name Model names
#' @param fig_path Figure output path
#' @param fig_name Figure name
#'
#' @return
#' @export
#'
#' @examples
compare_5sptl_models = function(model_path, model_name, fig_path, fig_name) {
  
pdf(here(fig_path, fig_name), width = 15)

require(here); require(tidyverse); require(expm); require(MASS); library(here)

# Read in models
model_list = list() # list to store models in 
model_sd = list() # list to store models in (sdrep)
for(i in 1:length(model_path)) model_list[[i]] = readRDS(file.path(paste(model_path[[i]], "/mle_report.RDS", sep = "")))
for(i in 1:length(model_path)) model_sd[[i]] = readRDS(file.path(paste(model_path[[i]], "/sd_report.RDS", sep = "")))

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
design_survey_index = design_survey_index %>% group_by(Country, Year, area_lab) %>% 
  summarise(sum_estimates = sum(area_RPN, na.rm = T), sum_var = sum(var_area_RPN, na.rm = T), 
            se = log_sigma(sqrt(sum_var)/sum_estimates), LCI = lognormal_CI(sum_estimates, se, 0.95)$lower, 
            UCI = lognormal_CI(sum_estimates, se, 0.95)$upper)



# Fits to tag data --------------------------------------------------------
# Compare fits to tag data
if(length(dim(model_list[[1]]$obs_tag_recovery)) == 3) {
    obs_tag_recovery = reshape2::melt(model_list[[1]]$obs_tag_recovery) %>% 
      rename(release_event = Var1, region = Var2, year = Var3, obs = value) %>%  # get observed
      mutate(age = 1e5) # arbitrary age set for left joining
}
  
if(length(dim(model_list[[1]]$obs_tag_recovery)) == 4) {
  obs_tag_recovery = reshape2::melt(model_list[[1]]$obs_tag_recovery) %>% 
    rename(age = Var1, release_event = Var2, region = Var3, year = Var4, obs = value) # get observed
}
  
# Get predicted tag recovery
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
  
  # Plot fits to data
print(  tag_recovery_df %>% 
          group_by(year, model, region) %>% 
          summarize(obs = sum(obs, na.rm = TRUE), pred = sum(pred, na.rm = TRUE)) %>% 
          ggplot() +
          geom_point(aes(x = year + 1977, y = obs), size = 2) +
          geom_line(aes(x = year + 1977, y = pred, color = model, lty = model), size = 1) +
          facet_wrap(~region, scales = "free", nrow = 3) +
          theme_bw() +
          labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Tag Recovery"))

# Get nLL
print(
  tag_recovery_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Tag Recovery nLL") +
    theme_bw()
)
  
## Fits to index data ------------------------------------------------------
# Get predicted tag recovery
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
         model = factor(model, levels = model_name))

# Plot survey biomass index residuals
print(
  ggplot(srv_bio_df) +
    geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
    geom_point(aes(x = Year, y = resids, color = model), size = 2, alpha = 0.75) +
    geom_smooth(aes(x = Year, y = resids, color = model, lty = model), se = FALSE, lwd = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Index Residuals") +
    theme_bw()
)

print(
  ggplot(srv_bio_df) +
    geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI)) +
    geom_line(aes(x = Year, y = pred, color = model),  lwd = 1) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Survey Fits") +
    theme_bw()
)

# Get nLL
print(
  srv_bio_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Index nLL") +
    theme_bw()
)

## Fits to comp data -------------------------------------------------------
### Fixed Gear (Age) -------------------------------------------------------------
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
    left_join(pred_fixed_age, by = c("Sex_Age", "region", "Year")) %>% 
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
  
  # plot fixed gear age average fits (region)
print(
  ggplot() +
    geom_col(fixed_age_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(region, Sex_Age, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.5) +
    geom_line(fixed_age_df %>% 
                group_by(region, Sex_Age, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    facet_wrap(~region, scales = "free") +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Age") +
    theme_bw()
)
  
  # plot fixed gear age average fits (all aggregated)
print(
  ggplot() +
    geom_col(fixed_age_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(Sex_Age, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.5) +
    geom_line(fixed_age_df %>% 
                group_by(Sex_Age, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Age") +
    theme_bw()
)

# Combined fixed gear age nLL
print(
  fixed_age_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Fixed Gear Age nLL") +
    theme_bw()
)
  
### Fixed Gear (Length) -------------------------------------------------------------
obs_fixed_len = reshape2::melt(model_list[[1]]$obs_fixed_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
    group_by(region, Year) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)
  
# get predicted fixed gear
pred_fixed_len = data.frame()
for(i in 1:length(model_list)) {
  pred_fixed_len_tmp = reshape2::melt(model_list[[i]]$pred_fixed_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[3])
  pred_fixed_len = rbind(pred_fixed_len, pred_fixed_len_tmp)
} # end i


fixed_len_df = obs_fixed_len %>% 
    left_join(pred_fixed_len, by = c("Sex_Len", "region", "Year")) %>% 
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
  
# plot fixed gear length average fits (region)
print(
  ggplot() +
    geom_col(fixed_len_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(region, Sex_Len, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.5) +
    geom_line(fixed_len_df %>% 
                group_by(region, Sex_Len, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    facet_wrap(~region, scales = "free") +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Length") +
    theme_bw()
)
  
# plot fixed gear length average fits (all aggregated)
print(
  ggplot() +
    geom_col(fixed_len_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(Sex_Len, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.5) +
    geom_line(fixed_len_df %>% 
                group_by(Sex_Len, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Length") +
    theme_bw()
)

# Combined fixed gear length nLL
print(
  fixed_len_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Fixed Gear Length nLL") +
    theme_bw()
)
  
### Trawl Gear (Length) -------------------------------------------------------------
obs_trwl_len = reshape2::melt(model_list[[1]]$obs_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
    group_by(region, Year) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)

# get predicted trawl gear
pred_trwl_len = data.frame()
for(i in 1:length(model_list)) {
  pred_trwl_len_tmp = reshape2::melt(model_list[[i]]$pred_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, pred = value) %>%  # get predicted 
    mutate(model = model_name[i], nll = model_list[[i]]$nll[2])
  pred_trwl_len = rbind(pred_trwl_len, pred_trwl_len_tmp)
} # end i
  
trwl_len_df = obs_trwl_len %>% 
    left_join(pred_trwl_len, by = c("Sex_Len", "region", "Year")) %>% 
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
  
  # plot trawl gear length average fits (region)
print(
  ggplot() +
    geom_col(trwl_len_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(region, Sex_Len, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.5) +
    geom_line(trwl_len_df %>% 
                group_by(region, Sex_Len, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    facet_wrap(~region, scales = "free") +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Trawl Gear Length") +
    theme_bw()
)
  
  # plot trawl gear length average fits (all aggregated)
print(
  ggplot() +
    geom_col(trwl_len_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(Sex_Len, model) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Len, y = Obs_Mean), alpha = 0.5) +
    geom_line(trwl_len_df %>% 
                group_by(Sex_Len, model) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Trawl Gear Length") +
    theme_bw()
)

# Combined trawl gear length nLL
print(
  trwl_len_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Trawl Gear Length nLL") +
    theme_bw()
)
  
### Survey Gear (age) -------------------------------------------------------------
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
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
      Country = ifelse(Country == 1, "Japan", "United States"),
      model = factor(model, levels = model_name)
    ) 
  
  
  # plot srv age length average fits (region)
print(
  ggplot() +
    geom_col(srv_age_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(region, Sex_Age, model, Country) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.5) +
    geom_line(srv_age_df %>% 
                group_by(region, Sex_Age, model, Country) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    facet_grid(region~Country, scales = "free") +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Srv Age") +
    theme_bw()
)
  
# plot srv gear age average fits (all aggregated)
print(
  ggplot() +
    geom_col(srv_age_df %>% 
               filter(model == model_name[1]) %>% 
               group_by(Sex_Age, model, Country) %>%
               summarize(Obs_Mean = mean(obs_prop)),
             mapping = aes(x = Sex_Age, y = Obs_Mean), alpha = 0.5) +
    geom_line(srv_age_df %>% 
                group_by(Sex_Age, model, Country) %>%
                summarize(Pred_Mean = mean(pred)),
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1) +
    facet_wrap(~Country, scales = "free") +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Srv Age") +
    theme_bw()
)

# Combined trawl gear length nLL
print(
  srv_age_df %>% 
    group_by(model) %>% 
    summarize(nll = mean(nll)) %>% 
    ggplot(aes(x = model, y = nll)) +
    geom_point(size = 5) +
    geom_line(group = 1, lty = 2, lwd = 1.2) +
    labs(x = "Model", y = "Combined Srv Age nLL") +
    theme_bw()
)

# Biomass and Recruitment -------------------------------------------------
# Read in models
ts_dat = data.frame()
for(i in 1:length(model_list)) {
  
  # get ssb
  ssb_tmp = reshape2::melt(model_list[[i]]$SSB_yr) %>% 
    rename(Year = Var1, region = Var2) %>% 
    mutate(Model = model_name[i], Type = "SSB")
  
  # get rec
  rec_tmp = reshape2::melt(model_list[[i]]$recruitment_yr) %>% 
    rename(Year = Var1, region = Var2) %>% 
    mutate(Model = model_name[i], Type = "Rec")

  # Get totals
  ssb_total_tmp = ssb_tmp %>% dplyr::select(-region) %>% 
    group_by(Year, Model, Type) %>% 
    summarize(value = sum(value)) %>% 
    mutate(region = 1e10)
  
  # Get totals
  rec_total_tmp = rec_tmp %>% dplyr::select(-region) %>% 
    group_by(Year, Model, Type) %>% 
    summarize(value = sum(value)) %>% 
    mutate(region = 1e10)

  ts_dat = rbind(ts_dat, ssb_tmp, rec_tmp, ssb_total_tmp, rec_total_tmp)
}

# Bind together
quants_df = ts_dat %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA", 
      region == 1e10 ~ "Total"
    ), 
    region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA", "Total")),
    Model = factor(Model, levels = model_name)
  )
    
print(
  ggplot(quants_df %>% filter(Type == "SSB"), aes(x = Year + 1959, y = value, color = Model)) +
    geom_line(lwd = 0.85) +
    facet_wrap(~region, scales = "free") +
    ylim(0, NA) +
    labs(x = "Year", y = "SSB", color = "Model", lty = "Model") +
    theme_bw()
)

print(
  ggplot(quants_df %>% filter(Type == "Rec")) +
    geom_line(aes(x = Year + 1959, y = value, color = Model), lwd = 0.85) +
    facet_wrap(~region) +
    labs(x = "Year", y = "SSB", color = "Model", lty = "Model") +
    theme_bw()
)

print(
  ggplot(quants_df %>% filter(Type == "Rec")) +
    geom_col(aes(x = Year + 1959, y = value, fill = Model), lwd = 0.85, position = "dodge") +
    facet_wrap(~region) +
    labs(x = "Year", y = "SSB", color = "Model", lty = "Model") +
    theme_bw()
)


# Movement ----------------------------------------------------------------
# get movement estimates
movement_model_est = data.frame()
stationary_movement_all = data.frame()
for(i in 1:length(model_list)) {
  
  # Loop through to get estimates and stationary distribution
  if(length(dim(model_list[[i]]$movement_matrix)) == 3) {
    
    stationary_dist = data.frame()
    for(t in 1:dim(model_list[[i]]$movement_matrix)[3]) {
      stationary_dist_mat = model_list[[i]]$movement_matrix[,,t] %^% 50 # get stationary distribution
      stationary_dist_tmp = data.frame(region = c("BS", "AI", "WGOA", "CGOA", "EGOA"), # munge into dataframe
                                       move = stationary_dist_mat[1,], Model = model_name[i], 
                                       TimeBlock = paste("TimeBlk", t), AgeBlock = paste("AgeBlk", 1))
      stationary_dist = rbind(stationary_dist, stationary_dist_tmp)
    } # end t
    
    # Extract and reformat estimates
    movement_est_tmp = reshape2::melt(model_list[[i]]$movement_matrix) %>% 
      rename(From = Var1, To = Var2, Time = Var3) %>% 
      mutate(Model = model_name[i], Age = paste("AgeBlk", 1), Time = paste("TimeBlk", Time))
  } # no age blocks
  
  if(length(dim(model_list[[i]]$movement_matrix)) == 4) {
    
    stationary_dist = data.frame()
    for(t in 1:dim(model_list[[i]]$movement_matrix)[3]) {
      for(a in 1:dim(model_list[[i]]$movement_matrix)[4]) {
        stationary_dist_mat = model_list[[i]]$movement_matrix[,,t,a] %^% 50 # get stationary distribution
        stationary_dist_tmp = data.frame(region = c("BS", "AI", "WGOA", "CGOA", "EGOA"), # munge into dataframe
                                         move = stationary_dist_mat[1,], Model = model_name[i], 
                                         TimeBlock = paste("TimeBlk", t), AgeBlock = paste("AgeBlk", a))
        stationary_dist = rbind(stationary_dist, stationary_dist_tmp)
      }
    } # end t
    
    
    movement_est_tmp = reshape2::melt(model_list[[i]]$movement_matrix) %>% 
      rename(From = Var1, To = Var2, Time = Var3, Age = Var4) %>% 
      mutate(Model = model_name[i], Time = paste("TimeBlk", Time),
             Age = paste("AgeBlk", Age))
  } # if there are age blocks
  
  movement_model_est = rbind(movement_model_est, movement_est_tmp)
  stationary_movement_all = rbind(stationary_movement_all, stationary_dist)
}

# residual munging
movement_df = movement_model_est %>% 
  mutate(
    From = case_when(
      From == 1 ~ "BS",
      From == 2 ~ "AI",
      From == 3 ~ "WGOA",
      From == 4 ~ "CGOA",
      From == 5 ~ "EGOA"
    ), From = factor(From, levels = rev(c("BS", "AI", "WGOA", "CGOA", "EGOA"))),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), To = factor(To, levels = (c("BS", "AI", "WGOA", "CGOA", "EGOA"))),
    Model = factor(Model, levels = model_name)
  )

for(i in 1:length(unique(movement_df$Model))) { 
  movement_plot_df = movement_df %>% filter(Model == model_name[i]) # filter to model
  print(
    ggplot(movement_plot_df, aes(x = To, y = From, fill = value, label = round(value, 3))) +
      geom_tile(alpha = 0.5) + 
      geom_text() +
      scale_fill_viridis_c() +
      facet_grid(Age~Time) +
      theme_test() +
      labs(fill = "Movement Probability", title = model_name[i])
  )
}

stationary_movement_all = stationary_movement_all %>% 
  mutate(
    region =  factor(region, levels = (c("BS", "AI", "WGOA", "CGOA", "EGOA"))),
    Model = factor(Model, levels = model_name)
  )

print(ggplot(stationary_movement_all, aes(x = region , y = move, color = Model, group = Model)) +
        geom_line(lwd = 1.3) +
        facet_grid(TimeBlock~AgeBlock) +
        theme_bw() +
        labs(x = "Region", y = "Stationary Movement Probability"))


# Reporting Rate ----------------------------------------------------------

# get predicted reporting rate
rep_rate_df = data.frame()
for(i in 1:length(model_list)) {
  rep_rate_tmp = reshape2::melt(model_list[[i]]$tag_reporting_rate) %>% 
    rename(region = Var1, Year = Var2) %>% mutate(Year = Year + 1977, model = model_name[i])
  rep_rate_df = rbind(rep_rate_df, rep_rate_tmp)
} # end i

# do residual munging
rep_rate_df = rep_rate_df %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), 
    region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA", "Total")),
    model = factor(model, levels = model_name)
  )

print(ggplot(rep_rate_df, aes(x = Year, y = value, color = model)) +
        geom_line(lwd = 1.3) +
        facet_wrap(~region, scales = "free") +
        theme_bw() +
        labs(x = "Year", y = "Reporting Rate", color = "Model"))


# Mean Recruitment --------------------------------------------------------
rec_par_df = data.frame()
for(i in 1:length(model_sd)) {
  pars = model_sd[[i]]$par.fixed[names(model_sd[[i]]$par.fixed) == "ln_mean_rec"] # get mle
  sd = sqrt(diag(model_sd[[1]]$cov.fixed)[names(diag(model_sd[[1]]$cov.fixed)) == "ln_mean_rec"]) # get se
  rec_tmp = data.frame(pars = pars, sd = sd, lwr = pars - 1.96 * sd, upr = pars + 1.96 * sd, model = model_name[i],
                       region = factor(c("BS", "AI", "WGOA", "CGOA", "EGOA"), levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) # munge into dataframe
  rec_par_df = rbind(rec_par_df, rec_tmp)
} # end i
 
print(
  ggplot(rec_par_df, aes(x = region, y = pars, ymin = lwr, ymax = upr, color = model)) +
    geom_pointrange(position = position_dodge(width = 0.3)) +
    theme_bw() +
    labs(x = "Region", y = "log(Mean Recruitment)", color = "Model")
)
  

# Catchability -----------------------------------------------------

q_par_df = data.frame()
for(i in 1:length(model_sd)) {
  pars = model_sd[[i]]$par.fixed[names(model_sd[[i]]$par.fixed) == "trans_srv_q"] # get mle
  sd = sqrt(diag(model_sd[[1]]$cov.fixed)[names(diag(model_sd[[1]]$cov.fixed)) == "trans_srv_q"]) # get se
  if(length(pars) > 2) {
    q_tmp = data.frame(pars = pars, sd = sd, lwr = pars - 1.96 * sd, upr = pars + 1.96 * sd, model = model_name[i],
                       region = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
                       srv = rep(c("jpll", "domll"), each = 5)) # munge into dataframe
  } else{
    q_tmp = data.frame(pars = pars, sd = sd, lwr = pars - 1.96 * sd, upr = pars + 1.96 * sd, model = model_name[i],
                       region = "Total", srv = c("jpll", "domll")) # munge into dataframe
  }
  q_par_df = rbind(q_par_df, q_tmp)
} # end i

# Munge region name
q_par_df = q_par_df %>% 
  mutate(region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA", "Total")))

print(
  ggplot(q_par_df, aes(x = region, y = pars, ymin = lwr, ymax = upr, color = model)) +
    geom_pointrange(position = position_dodge(width = 0.3)) +
    theme_bw() +
    facet_wrap(~srv) +
    labs(x = "Region", y = "log(Survey Q)", color = "Model")
)
  

# Selectivity -------------------------------------------------------------



dev.off()
  
}






