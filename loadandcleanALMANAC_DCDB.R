library(synergyfinder)
library(data.table)
library(dplyr)
library(tidyverse)
library("GGally")
library(viridis)

#dcdb calculation 
data_dir = getwd()
dcdb_score <- read.csv(paste0(data_dir, "/dcdb_crawler_almanac_withID.csv"))
NCI60_TO_TISSUE <- read.csv(paste0(data_dir, "/NCI60_TO_TISSUE.txt"),sep = '\t')

joined_tibble <- left_join(dcdb_score, NCI60_TO_TISSUE, 
                      by = c("cellName" = "Cell.Line.Name"))

long_NCI60_COMBO_SCORE <- read.csv(paste0(data_dir, "/long_DTP_NCI60_ALMANAC_COMBO_SCORE_withID.csv"))
combosyn_5scores_temp1 <- left_join(joined_tibble, long_NCI60_COMBO_SCORE,
                              by = c("pubchemID_x" = "pubchemID_x", "pubchemID_y" = "pubchemID_y",
                                     "cellName"="CELLNAME"))
combosyn_5scores_temp1 <- combosyn_5scores_temp1[!is.na(combosyn_5scores_temp1$COMBO_SCORE),]
combosyn_5scores_temp2 <- left_join(joined_tibble, long_NCI60_COMBO_SCORE,
                                    by = c("pubchemID_y" = "pubchemID_x", "pubchemID_x" = "pubchemID_y",
                                           "cellName"="CELLNAME"))
combosyn_5scores_temp2 <- combosyn_5scores_temp2[!is.na(combosyn_5scores_temp2$COMBO_SCORE),]

combosyn_5scores <- rbind(combosyn_5scores_temp1,combosyn_5scores_temp2) %>% 
  distinct(pubchemID_x, pubchemID_y,cellName, .keep_all = T)
rm(combosyn_5scores_temp1)
rm(combosyn_5scores_temp2)

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
    geom_point(aes(colour=density,  alpha = 0.2, size=3), ...) +
    scale_color_viridis()      
  p
}
fig1c <- ggpairs(combosyn_5scores,
                 columns=c("Loewe","ZIP", "Bliss", "HSA","COMBO_SCORE"),
                 columnLabels=c("Loewe", "ZIP","Bliss", "HSA","ALMANAC.Score"),
                 title = "Correlation between synergy matrics",
                 lower=list(continuous=my_fn),
                 diag = list(continuous = "barDiag", fill="green",colour = "green"),
                 upper = list(continuous = wrap("cor", method = "pearson"))) +
  mytheme
fig1c


