#'
#' Compare removals from 5 area model between single F apportioned by survey
#' VS single area reference points 
#'

source("(00) Init.R")
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(SpatialSablefishAssessment)
library(TMB)


mod_5A = file.path(DIR$app_5A,model_labels, "Model1960_03_srva")

