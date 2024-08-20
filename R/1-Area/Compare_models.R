#'
#' Compare all 1 area models using bookdown code
#'
source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)
## run 1977 inital model developments
model_labels = c("1977_starta", "1977_init_devsa", "1977_init_num_devsa", "sum_zero_recruita", "DM_comp_likelihoodsa")#, "sum_zero_recruita", "include_tag_data", "negative_binomiala","decadal_tag_reporta")
bookdown_labels = c("1977 Finit","1977 Finit & 15 ninit", '1977 Finit & 28 ninit', "1977 Finit sum 0", "1977 DM")#, "tag", "NB", "decade-tag")
model_dir = file.path(DIR$app_1A,model_labels)
file.exists(model_dir)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1977start")
model_description = '
- "1977 Finit"  Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial.
- "1977 Finit & 15 ninit" Start year is 1977, estimate an initial F-init parameter and 15 initial devs for ages 3->18 with ages 18->31 sharing the last estimated value. Comp data assumed multinomial.
- "1977 Finit & 28 ninit" Start year is 1977, estimate an initial F-init parameter and 28 initial devs for ages 3->30. Comp data assumed multinomial.
- "1977 Finit sum 0" Start year 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial. Constrain Recruitment deviations to sum = 0 over the entire time-series.
- "1977 DM"  Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed to be Dirichlet-multinomial.
'
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_together_1977start")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

## run 1977 tag models
model_labels = c("1977_starta", "include_tag_dataa", "negative_binomiala", "decadal_tag_reporta")#, "sum_zero_recruita", "include_tag_data", "negative_binomiala","decadal_tag_reporta")
bookdown_labels = c("1977 Finit","1977 Poisson", 'NB', "NB TR decade")#, "tag", "NB", "decade-tag")
model_dir = file.path(DIR$app_1A,model_labels)
file.exists(model_dir)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1977Tag")
model_description = '
- "1977 Finit"  Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial.
- "1977 Poisson" Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial. Included tagging data with Poisson likelihood, estimated time-invariant tag-reporting rate, assumed constant across space
- "NB" Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial. Included tagging data with Negative binomial likelihood, estimated time-invariant tag-reporting rate, assumed constant across space
- "NB TR decade" Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial. Included tagging data with Negative binomial likelihood, estimated tag-reporting rate for each decade, assumed constant across space
'
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_together_1977Tag")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)


## compare 1960 models
model_labels = c("1960_starta", "1960_no_inita", "1960_dont_est_10a", "1960_sum_zero_recruita","1960_include_tag_dataa", "1960_negative_binomiala","1960_DM_comp_likelihoodsa" )#, "Model_02a")#, "1960_TagData_01a", "1960_TagData_02a","1960_TagData_03a")
bookdown_labels = c("1960", "no-init dev", "fixed 10 YCS","sum zero rec", "tag-pois", "tag-NB", "decade-TR")#, "NB reg rec")#, "Poisson tag","NB tag", "decade TR", "reduced_initage)
model_description = '
- "1960"  no movement, no tagging data, comp likelihood assumed multinomial, global rec devs, estimate F-init (not intial age-devs)
- "no-init dev"  dont estimate the initial age deviations
- "fixed 10 YCS" fix the first 10 rec-devs to be 0
- "sum zero rec" Constrain recruitment deviations to sum 0 using the simplex
- "tag-pois" include tagging data which is assumed to be Poisson
- "tag-NB" tagging assumed to be negative binomial with estimated overdispersion parameter
- "decade-TR" add decadal estimated tag-reporting parameters
'
model_dir = file.path(DIR$app_1A,model_labels)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1960_init")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_together_1960_init")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

## Compare across start year
## 1960 models
model_labels = c("decadal_tag_reporta", "1960_DM_comp_likelihoodsa")#, "1960_include_tag_dataa", "1960_negative_binomiala","1960decadal_tag_reporta")
model_dir = file.path(DIR$app_1A,model_labels)
file.exists(model_dir)
bookdown_labels = c( "1977", "1960")#, "tag_data","negative_binomiala","decadal_tag_reporta")
model_description = '
- "1977"  Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial.
- "1960" Start year is 1960, assumes age-structure was in equilibrium. Comp data assumed DM
'

bookdown_dir = file.path(DIR$book1A, "Comparison_together_1960_1977_init")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1960_1977_init")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

## 1960 models with Tag-data
model_labels = c("1977_init_devsa", "1960_sum_zero_recruita", "decadal_tag_reporta", "1960_DM_comp_likelihoodsa")#, "1960_include_tag_dataa", "1960_negative_binomiala","1960decadal_tag_reporta")
model_dir = file.path(DIR$app_1A,model_labels)
file.exists(model_dir)
bookdown_labels = c( "1977", "1960", "1977 tag", "1960 tag")#, "tag_data","negative_binomiala","decadal_tag_reporta")
model_description = '
- "1977"  Start year is 1977, estimate an initial F-init parameter assumed age-structure was in equilibrium. Comp data assumed multinomial.
- "1960" Start year is 1960, assumes age-structure was in equilibrium. Comp data assumed multinomial. Did not estimate the first 10 recrutiment deviations as a sensitivity
- "Poisson" Same as "1960", but included tag-recapture data with an assumed Poisson likelihood
- "NB" Same as "Poisson", but assumed tag-recapture likelihood was negative binomial
- "NB decade TR" Same as "NB", but estimated a tag-reporting parameter for each decade.
'

bookdown_dir = file.path(DIR$book1A, "Comparison_together_1977_1960_tag")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1977_1960_tag")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)

model_labels = c("1977_init_num_devsa", "1960_starta",  "1960_dont_est_10a", "1960_sum_zero_recruita")
model_dir = file.path(DIR$app_1A,model_labels)
bookdown_labels = c( "1977", "1960", "1960 est YCS", "sum_zero")
bookdown_dir = file.path(DIR$book1A, "Comparison_together_1960_init")
summarise_multiple_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
bookdown_dir = file.path(DIR$book1A, "Comparison_individual_1960_init")
summarise_individual_models(model_dir = model_dir, bookdown_labels = bookdown_labels, bookdown_dir = bookdown_dir, model_description = model_description)
