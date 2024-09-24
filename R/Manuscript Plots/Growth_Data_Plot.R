# Purpose: To visualize growth data
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date 9/23/24

# Set up ------------------------------------------------------------------
library(tidyverse)
library(here)

# Read in age data
age_dat = read.csv(here('Data',  'Survey', 'age_view.csv'), check.names = FALSE)

# Clean up some age data
age_dat = age_dat %>% 
filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
       `Sex Description` != "Unknown", `Error Flag` == 0) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`,
         Haul = `Station Number`,
         Mgmt_Area = `NPFMC Sablefish Mgmt Area`) %>% 
  filter(!Age %in% c(0, 1)) %>% 
  select(Sex_name, Length, Weight, Age, Year, Haul, Mgmt_Area) %>% 
  mutate(Sex_name = ifelse(Sex_name == "female", "Females", 'Males'),
         Mgmt_Area = ifelse(Mgmt_Area == "Aleutians", "Aleutian Islands", Mgmt_Area))


# Plot! -------------------------------------------------------------------

growth_plot = ggplot(age_dat, aes(x = Age, y = Weight, color = Mgmt_Area)) +
  geom_point(alpha = 0.05, size = 2.5) + 
  geom_smooth(method = 'loess', lwd = 1.3) +
  scale_color_manual(values = c("black", "#56B4E9", "#009E73","#0072B2", "#D55E00", "#E69F00")) +
  facet_wrap(~Sex_name) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top')  +
  labs(x = "Age", y = "Weight (kg)", color = "Sablefish Management Regions")

ggsave(
  here("Figs", "Manuscript_Plots", "Growth.png"),
  growth_plot, width = 15
)
