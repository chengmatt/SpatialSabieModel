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
compare_2_models = function(model_path, model_name, fig_path, fig_name) {
  
  pdf(here(fig_path, fig_name), width = 13)
  
  require(here); require(tidyverse)
  
  # Read in models
  model_1 = readRDS(file.path(paste(model_path[[1]], "/mle_report.RDS", sep = "")))
  model_2 = readRDS(file.path(paste(model_path[[2]], "/mle_report.RDS", sep = "")))
  
  # Read in abundance index used
  design_survey_index = readRDS(file = file.path("Data", "Survey", "regional_abundance_estimates.RDS"))
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
  obs_tag_recovery = reshape2::melt(model_1$obs_tag_recovery) %>% 
    rename(release_event = Var1, region = Var2, year = Var3, obs = value) # get observed
  pred_tag_recovery_01 = reshape2::melt(model_1$pred_tag_recovery) %>% 
    rename(release_event = Var1, region = Var2,  year = Var3, Mod1 = value) 
  pred_tag_recovery_02 = reshape2::melt(model_2$pred_tag_recovery) %>% 
    rename(release_event = Var1, region = Var2, year = Var3, Mod2 = value) 
  
  # Join all datasets together
  tag_recovery_df = obs_tag_recovery %>% 
    left_join(pred_tag_recovery_01, by = c("release_event", "region", "year")) %>% 
    left_join(pred_tag_recovery_02, by = c("release_event", "region", "year")) %>% 
    pivot_longer(cols = c("Mod1", "Mod2"), values_to = "pred", names_to = 'model') %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    ) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
  # Plot fits to data
print(  tag_recovery_df %>% 
          group_by(year, model, region) %>% 
          summarize(obs = sum(obs), pred = sum(pred)) %>% 
          ggplot() +
          geom_point(aes(x = year + 1977, y = obs), size = 2) +
          geom_line(aes(x = year + 1977, y = pred, color = model, lty = model), size = 1.3) +
          facet_wrap(~region, scales = "free", nrow = 3) +
          theme_bw() +
          labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Tag Recovery"))
  
## Fits to index data ------------------------------------------------------
  pred_srv_01 = reshape2::melt(model_1$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, Mod1 = value) %>%  
    mutate(Year = Year + 1959, 
           area_lab = case_when(
             area_lab == 1 ~ "BS",
             area_lab == 2 ~ "AI",
             area_lab == 3 ~ "WGOA",
             area_lab == 4 ~ "CGOA",
             area_lab == 5 ~ "EGOA"),
           Country = ifelse(Country == 1, "Japan", "United States")) %>% 
    filter(Year >= 1979)
  
  pred_srv_02 = reshape2::melt(model_2$pred_srv_bio) %>% 
    rename(area_lab = Var1, Year = Var2, Country = Var3, Mod2 = value) %>%  
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
    left_join(cbind(pred_srv_01,  Mod2 = pred_srv_02$Mod2), by = c("Country", "area_lab", "Year")) %>% 
    pivot_longer(cols = c('Mod1', "Mod2"), values_to = "Pred", names_to = "model") %>% 
    mutate(resids = (log(sum_estimates) - log(Pred)) / se) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
  # Plot survey biomass
print(
  ggplot(srv_bio_df) +
    geom_segment(aes(x = Year, y = 0, xend = Year, yend = resids), lwd = 0.2) +
    geom_point(aes(x = Year, y = resids, color = model), size = 3, alpha = 0.85) +
    geom_smooth(aes(x = Year, y = resids, color = model, lty=  model), se = FALSE, lwd = 1.3) +
    geom_hline(yintercept = 0, lty = 2) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Index Residuals") +
    theme_bw()
)

print(
  ggplot(srv_bio_df) +
    geom_pointrange(aes(x = Year, y = sum_estimates, ymin = LCI, ymax = UCI)) +
    geom_line(aes(x = Year, y = Pred, color = model, lty=  model),  lwd = 1.3) +
    facet_grid(area_lab~Country, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Survey Fits") +
    theme_bw()
)
  
  ## Fits to comp data -------------------------------------------------------
  ### Fixed Gear (Age) -------------------------------------------------------------
  obs_fixed_age = reshape2::melt(model_1$obs_fixed_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value) %>% 
    group_by(region, Year) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)
  
  pred_fixed_age_01 = reshape2::melt(model_1$pred_fixed_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Mod1 = value) # get predicted 
  pred_fixed_age_02 = reshape2::melt(model_2$pred_fixed_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Mod2 = value) # get predicted 
  
  fixed_age_df = obs_fixed_age %>% 
    left_join(cbind(pred_fixed_age_01, Mod2 = pred_fixed_age_02$Mod2),
              by = c("Sex_Age", "region", "Year")) %>% 
    pivot_longer(cols = c("Mod1", "Mod2"), values_to = "pred", names_to = 'model') %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    ) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
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
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
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
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Age") +
    theme_bw()
)
  
  ### Fixed Gear (Length) -------------------------------------------------------------
  obs_fixed_len = reshape2::melt(model_1$obs_fixed_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
    group_by(region, Year) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)
  
  pred_fixed_len_01 = reshape2::melt(model_1$pred_fixed_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Mod1 = value) # get predicted 
  pred_fixed_len_02 = reshape2::melt(model_2$pred_fixed_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Mod2 = value) # get predicted 
  
  fixed_len_df = obs_fixed_len %>% 
    left_join(cbind(pred_fixed_len_01, Mod2 = pred_fixed_len_02$Mod2),
              by = c("Sex_Len", "region", "Year")) %>% 
    pivot_longer(cols = c("Mod2", "Mod1"), values_to = "pred", names_to = 'model') %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    ) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
  
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
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
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
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Fixed Gear Length") +
    theme_bw()
)
  
  ### Trawl Gear (Length) -------------------------------------------------------------
  obs_trwl_len = reshape2::melt(model_1$obs_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Obs = value) %>% 
    group_by(region, Year) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)
  
  pred_trwl_len_01 = reshape2::melt(model_1$pred_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Mod1 = value) # get predicted 
  pred_trwl_len_02 = reshape2::melt(model_2$pred_trwl_catchatlgth) %>% 
    rename(Sex_Len = Var1, region = Var2, Year = Var3, Mod2 = value) # get predicted 
  
  trwl_len_df = obs_trwl_len %>% 
    left_join(cbind(pred_trwl_len_01, Mod2 = pred_trwl_len_02$Mod2),
              by = c("Sex_Len", "region", "Year")) %>% 
    pivot_longer(cols = c("Mod2", "Mod1"), values_to = "pred", names_to = 'model') %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    ) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
  
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
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
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
              mapping = aes(x = Sex_Len, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
    labs(x = "Sex-Length Category", y = "Value", color = "Model", lty = "Model", title = "Trawl Gear Length") +
    theme_bw()
)
  
  ### Survey Gear (age) -------------------------------------------------------------
  obs_srv_age = reshape2::melt(model_1$obs_srv_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Obs = value, Country = Var4) %>% 
    group_by(region, Year, Country) %>% 
    mutate(total = sum(Obs), obs_prop = Obs / total) %>% 
    filter(total != 0)
  
  srv_age_len_01 = reshape2::melt(model_1$pred_srv_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Mod1 = value, Country = Var4) 
  srv_age_len_02 = reshape2::melt(model_2$pred_srv_catchatage) %>% 
    rename(Sex_Age = Var1, region = Var2, Year = Var3, Mod2 = value, Country = Var4) 
  
  srv_age_df = obs_srv_age %>% 
    left_join(cbind(srv_age_len_01, Mod2 = srv_age_len_02$Mod2),
              by = c("Sex_Age", "region", "Year", "Country")) %>% 
    pivot_longer(cols = c("Mod2", "Mod1"), values_to = "pred", names_to = 'model') %>% 
    mutate(
      region = case_when(
        region == 1 ~ "BS",
        region == 2 ~ "AI",
        region == 3 ~ "WGOA",
        region == 4 ~ "CGOA",
        region == 5 ~ "EGOA"
      ), region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
      Country = ifelse(Country == 1, "Japan", "United States")
    ) %>% 
    mutate(model = ifelse(model == "Mod1", model_name[1], model_name[2]))
  
  
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
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
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
              mapping = aes(x = Sex_Age, y = Pred_Mean, color = model, lty = model), lwd = 1.5) +
    facet_wrap(~Country, scales = "free") +
    labs(x = "Sex-Age Category", y = "Value", color = "Model", lty = "Model", title = "Srv Age") +
    theme_bw()
  
)

# Biomass and Recruitment -------------------------------------------------
# Read in models
model_1_sd = readRDS(file.path(paste(model_path[[1]], "/sd_report.RDS", sep = "")))
model_2_sd = readRDS(file.path(paste(model_path[[2]], "/sd_report.RDS", sep = "")))

# From model 1
ssb_01 = reshape2::melt(matrix(model_1_sd$value[names(model_1_sd$value) %in% c("SSB_yr")], nrow = length(model_1$years), ncol = model_1$n_regions)) %>%  # get ssb
  mutate(Model = model_name[1], Type = "SSB") %>% rename(Year = Var1, region = Var2)
rec_01 = reshape2::melt(matrix(model_1_sd$value[names(model_1_sd$value) %in% c("recruitment_yr")], nrow = length(model_1$years), ncol = model_1$n_regions)) %>%  # get recruitment
  mutate(Model = model_name[1], Type = "Recruitment") %>%  rename(Year = Var1, region = Var2)

# From model 2
ssb_02 = reshape2::melt(matrix(model_2_sd$value[names(model_2_sd$value) %in% c("SSB_yr")], nrow = length(model_2$years), ncol = model_2$n_regions)) %>%  # get ssb
  mutate(Model = model_name[2], Type = "SSB") %>% rename(Year = Var1, region = Var2,)
rec_02 = reshape2::melt(matrix(model_2_sd$value[names(model_2_sd$value) %in% c("recruitment_yr")], 
                               nrow = length(model_2$years), ncol = model_2$n_regions)) %>%  # get recruitment
  mutate(Model = model_name[2], Type = "Recruitment") %>% rename(Year = Var1, region = Var2)

# Bind together
quants_df = rbind(ssb_01, rec_01, 
                  ssb_02, rec_02) %>% 
  mutate(
    region = case_when(
      region == 1 ~ "BS",
      region == 2 ~ "AI",
      region == 3 ~ "WGOA",
      region == 4 ~ "CGOA",
      region == 5 ~ "EGOA"
    ), 
    region = factor(region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))
    
print(
  ggplot(quants_df, aes(x = Year + 1959, y = value, color = Model, lty = Model)) +
    geom_line(lwd = 1.5) +
    facet_grid(Type~region, scales = "free") +
    labs(x = "Year", y = "Value", color = "Model", lty = "Model", title = "Time Series") +
    theme_bw()
)

dev.off()
  
}






