library(synergyfinder)
library(data.table)
library(dplyr)
library(tidyverse)
library("GGally")
library(viridis)

data_dir = getwd()
# Load the ALMANAC data.
almanac_dat = read.csv(paste0(data_dir, "/NCI-ALMANAC_full_data.csv"))
# Convert it to a data.table. 
setDT(almanac_dat)

# create a block_id for groupby
almanac_dat <- almanac_dat %>% group_by(Drug1,Drug2,CellLine) %>% mutate(block_id = cur_group_id())
screening_data <- almanac_dat %>% 
  dplyr::select(block_id,Drug1,Drug2,Conc1, Conc2, PercentageGrowth) %>% 
  rename( "drug_row" = "Drug1", 
          "drug_col" = "Drug2",
          "conc_r" = "Conc1",
          "conc_c" = "Conc2", 
          "response"="PercentageGrowth")

screening_data$conc_r <- screening_data$conc_r * 10^9
screening_data$conc_c <- screening_data$conc_c * 10^9
screening_data$conc_r_unit <- "nM"
screening_data$conc_c_unit <- "nM"
rm(almanac_dat)

#循环
empty_df = data.frame()

calsyn <- function(st,ed){
    test1 <- screening_data %>% filter(block_id %in% c(seq(st,ed))) #321240
    res <- synergyfinder::ReshapeData(
      data = test1,
      data_type = "viability",
      impute = TRUE,
      impute_method = NULL
    )
    res <- synergyfinder::CalculateSynergy(
      data = res,
      method = c("ZIP", "HSA", "Bliss", "Loewe"),
      Emin = NA,
      Emax = NA,
      correct_baseline = "non")
    
    empty_df <- rbind(empty_df,res$drug_pairs)
    empty_df
}

# >>>> Excute to calculate synergy score for each block of dose-response matrics
# d <- 321001:321240 #1:321240
# lgth <- 1 # split length
# chunks <- split(d, ceiling(seq_along(d)/lgth))
# 
# for (ck in chunks){
#   st = head(ck,1)
#   end = tail(ck,1)
#   
#   tryCatch(
#     expr = {
#       print(st)
#       empty_df <- calsyn(st,end)
#       message("Successfully executed the log(x) call.")
#       return(empty_df)
#     },
#     error = function(e){
#       message('Caught an error!')
#       print(e)
#     })
# }
# write.table(empty_df, file = "ComboDrugSyn_Nov2022.csv",sep=",")

combosyn <- read.csv(paste0(data_dir, "/ComboDrugSyn_Nov2022.csv"))
NCI60_TO_TISSUE <- read.csv(paste0(data_dir, "/NCI60_TO_TISSUE.txt"),sep = '\t')

# fig1b, map blockid to cell/tissue
map2cell <- screening_data %>% 
  dplyr::select(CellLine,block_id) %>% distinct(block_id, .keep_all = TRUE) 

joined_tibble <- left_join(combosyn, map2cell, 
                           by = c("block_id" = "block_id"))
combosyn <- left_join(joined_tibble, NCI60_TO_TISSUE, 
                           by = c("CellLine" = "Cell.Line.Name"))
#get almanac combo score
NCI60_COMBO_SCORE <- read.csv(paste0(data_dir, "/DTP_NCI60_ALMANAC_COMBO_SCORE.csv"),check.names=FALSE)
long_NCI60_COMBO_SCORE <- NCI60_COMBO_SCORE %>% pivot_longer(cols = -c(NSC1, NSC2), names_to = "CELLNAME", values_to = "COMBO_SCORE")

mylookup_tibble = read.csv(paste0(data_dir, "/ComboCompoundNames_small.txt"),sep='\t',)
joined_tibble <- left_join(long_NCI60_COMBO_SCORE, mylookup_tibble, 
                           by = c("NSC1" = "IX"))
long_NCI60_COMBO_SCORE <- left_join(joined_tibble, mylookup_tibble, 
                         by = c("NSC2" = "IX"))
write.table(long_NCI60_COMBO_SCORE, file = "long_DTP_NCI60_ALMANAC_COMBO_SCORE.csv",sep=",")

# >>>> finally get complete scores

#dcdb calculation 
dcdb_score <- read.csv(paste0(data_dir, "/dcdb_crawler_almanac_withID.csv"))
NCI60_TO_TISSUE <- mutate_all(NCI60_TO_TISSUE, .funs=toupper)
dcdb_score <- mutate_all(dcdb_score, .funs=toupper)
long_NCI60_COMBO_SCORE <- mutate_all(long_NCI60_COMBO_SCORE, .funs=toupper)

dcdb_score <- left_join(dcdb_score, NCI60_TO_TISSUE, 
                        by = c("Cell.line" = "Cell.Line.Name"))

combosyn_5scores <- left_join(dcdb_score, long_NCI60_COMBO_SCORE,
                                    by = c("Drug1" = "Drug.x", "Drug2" = "Drug.y",
                                           "Cell.line"="CELLNAME")) %>%
                              dplyr::select(ID,Drug1,Drug2,Cell.line,Panel.Name,
                                          ZIP,HSA,Bliss,Loewe,COMBO_SCORE)
# combosyn_5scores <- left_join(combosyn, dcdb_score,
#                                     by = c("drug1" = "Drug.x", "drug2" = "Drug.y",
#                                            "CellLine"="CELLNAME")) %>%
#                               dplyr::select(block_id,drug1,drug2,CellLine,Panel.Name,
#                                           ZIP_synergy,HSA_synergy,Bliss_synergy,Loewe_synergy,COMBO_SCORE)

# combosyn_5scores <- left_join(combosyn, long_NCI60_COMBO_SCORE, 
#                                     by = c("drug1" = "Drug.x", "drug2" = "Drug.y",
#                                            "CellLine"="CELLNAME")) %>% 
#                               dplyr::select(block_id,drug1,drug2,CellLine,Panel.Name,
#                                           ZIP_synergy,HSA_synergy,Bliss_synergy,Loewe_synergy,COMBO_SCORE)
combosyn_5scores <- combosyn_5scores[!is.na(combosyn_5scores$COMBO_SCORE),] # remove combo score is 0

combosyn_5scores$Loewe <- as.double(combosyn_5scores$Loewe)
combosyn_5scores$ZIP <- as.double(combosyn_5scores$ZIP)
combosyn_5scores$Bliss <- as.double(combosyn_5scores$Bliss)
combosyn_5scores$HSA <- as.double(combosyn_5scores$HSA)
combosyn_5scores$COMBO_SCORE <- as.double(combosyn_5scores$COMBO_SCORE)
fig1c <- ggpairs(combosyn_5scores,
                 columns=c("Loewe","ZIP", "Bliss", "HSA","COMBO_SCORE"),
                 columnLabels=c("Loewe", "ZIP","Bliss", "HSA","ALMANAC.Score"),
                 title = "Scatterplot Matrix with Linear and Loess Fits",
                 lower=list(continuous=my_fn),
                 diag = list(continuous = "barDiag", fill="green",colour = "green"),
                 upper = list(continuous = wrap("cor", method = "pearson"))) +
  mytheme
fig1c

#fig1c

mytheme <-  theme(strip.background = element_blank(),
                  panel.grid       = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(color="grey20", fill=NA))

my_fn <- function(data, mapping, N=100, ...){
  
  get_density <- function(x, y, n ) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  X <- eval_data_col(data, mapping$x)
  Y <- eval_data_col(data, mapping$y)
  
  data$density <- get_density(x=X, y=Y, n=N)
  
  p <- ggplot(data, mapping) +
    geom_point(aes(colour=density,  alpha = 0.2, size=5), ...) +
    scale_color_viridis()      
  p
}

#filter combosyn_5scores: Loewe, ZIP, Bliss, HSA on
combosyn_5scores <- combosyn_5scores %>% filter(Loewe_synergy > -200 & Loewe_synergy < 200)

fig1c <- ggpairs(combosyn_5scores,
        columns=c("Loewe_synergy","ZIP_synergy", "Bliss_synergy", "HSA_synergy","COMBO_SCORE"),
        columnLabels=c("Loewe", "ZIP","Bliss", "HSA","ALMANAC.Score"),
        title = "Scatterplot Matrix with Linear and Loess Fits",
        lower=list(continuous=my_fn),
        diag = list(continuous = "barDiag", fill="green",colour = "green"),
        upper = list(continuous = wrap("cor", method = "pearson"))) +
  mytheme
fig1c

#fig1d








