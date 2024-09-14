#' Title Compare 2 models
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
compare_sptl_models = function(model_path, model_name, fig_path, fig_name) {
  
  pdf(here(fig_path, fig_name), width = 15)
  
  require(here); require(tidyverse); require(expm); require(MASS); library(here)
  
  # Read in models
  model_list = list() # list to store models in 
  for(i in 1:length(model_path)) model_list[[i]] = readRDS(file.path(paste(model_path[[i]], "/mle_report.RDS", sep = "")))
  
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
  obs_tag_recovery = reshape2::melt(model_list[[1]]$obs_tag_recovery) %>% 
    rename(release_event = Var1, region = Var2, year = Var3, obs = value) # get observed
  
  # Get predicted tag recovery
  pred_tag_recovery = data.frame()
  for(i in 1:length(model_list)) {
    pred_tag_recovery_tmp = reshape2::melt(model_list[[i]]$pred_tag_recovery) %>% 
      rename(release_event = Var1, region = Var2,  year = Var3, pred = value) %>% 
      mutate(model = model_name[i])
    pred_tag_recovery = rbind(pred_tag_recovery, pred_tag_recovery_tmp)
  } # end i

  # Join all datasets together
  tag_recovery_df = obs_tag_recovery %>% 
    left_join(pred_tag_recovery, by = c("release_event", "region", "year")) %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )
  
  # Plot fits to data
print(  tag_recovery_df %>% 
          group_by(year, model, region) %>% 
          summarize(obs = sum(obs), pred = sum(pred)) %>% 
          ggplot() +
          geom_point(aes(x = year + 1977, y = obs), size = 2) +
          geom_line(aes(x = year + 1977, y = pred, color = model, lty = model), size = 1) +
          facet_wrap(~region, scales = "free", nrow = 3) +
          theme_bw() +
          labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Tag Recovery"))
  
## Fits to index data ------------------------------------------------------
# Get predicted tag recovery
pred_srv = data.frame()
for(i in 1:length(model_list)) {
  pred_srv_tmp = reshape2::melt(model_list[[i]]$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, pred = value) %>%  
    mutate(model = model_name[i])
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
  mutate(resids = (log(sum_estimates) - log(pred)) / se)

# Plot survey biomass
print(
  ggplot(srv_bio_df) +
    geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
    geom_point(aes(x = Year, y = resids, color = model), size = 1.5, alpha = 0.75) +
    geom_smooth(aes(x = Year, y = resids, color = model, lty=  model), se = FALSE, lwd = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Index Residuals") +
    theme_bw()
)

print(
  ggplot(srv_bio_df) +
    geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI)) +
    geom_line(aes(x = Year, y = pred, color = model, lty=  model),  lwd = 1) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Survey Fits") +
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
    mutate(model = model_name[i])
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
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
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
    mutate(model = model_name[i])
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
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
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
    mutate(model = model_name[i])
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
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
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
    mutate(model = model_name[i])
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
      Country = ifelse(Country == 1, "Japan", "United States")
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
    region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA", "Total")))
    
print(
  ggplot(quants_df %>% filter(Type == "SSB"), aes(x = Year + 1959, y = value, color = Model, lty = Model)) +
    geom_line(lwd = 1) +
    facet_wrap(~region, scales = "free") +
    labs(x = "Year", y = "SSB", color = "Model", lty = "Model", title = "Time Series") +
    theme_bw()
)

print(
  ggplot(quants_df %>% filter(Type == "Rec"), aes(x = Year + 1959, y = value, fill = Model)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~region, scales = "free") +
    labs(x = "Year", y = "SSB", color = "Model", lty = "Model", title = "Time Series") +
    theme_bw()
)

# Movement ----------------------------------------------------------------
# get movement estimates
movement_model_est = data.frame()
for(i in 1:length(model_list)) {
  movement_est_tmp = reshape2::melt(model_list[[i]]$movement_matrix) %>% 
    rename(From = Var1, To = Var2, Block = Var3) %>% 
    mutate(Model = model_name[i])
  movement_model_est = rbind(movement_model_est, movement_est_tmp)
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
    ), From = factor(From, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
    To = case_when(
      To == 1 ~ "BS",
      To == 2 ~ "AI",
      To == 3 ~ "WGOA",
      To == 4 ~ "CGOA",
      To == 5 ~ "EGOA"
    ), To = factor(To, levels = rev(c("BS", "AI", "WGOA", "CGOA", "EGOA"))))

print(
  ggplot(movement_df, aes(x = From, y = To, fill = value, label = round(value, 2))) +
    geom_tile(alpha = 0.5) + 
    geom_text() +
    scale_fill_viridis_c() +
    facet_grid(Block~Model) +
    theme_test() +
    labs(fill = "Movement Probability")
)

dev.off()
  
}






