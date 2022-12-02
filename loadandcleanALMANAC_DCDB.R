library(synergyfinder)
library(data.table)
library(dplyr)
library(tidyverse)
library("GGally")
library(viridis)

#dcdb calculation 
dcdb_score <- read.csv(paste0(data_dir, "/dcdb_crawler_almanac_withID.csv"))
NCI60_TO_TISSUE <- read.csv(paste0(data_dir, "/NCI60_TO_TISSUE.txt"),sep = '\t')

joined_tibble <- left_join(dcdb_score, NCI60_TO_TISSUE, 
                      by = c("cellName" = "Cell.Line.Name"))

