library(synergyfinder)
library(data.table)
library(dplyr)
library(tidyverse)
library("GGally")
library(viridis)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)

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
tiff("plot.png",units="in",width=5,height=5,res=300)
print(fig1c)
dev.off


#“synergistic” for each metric if their scores were within the top 5% and “antagonistic” if the scores were in the bottom 66.67%.

#row names for synergistic
assigning <- function(name) {
  print(name)
  print(quantile(combosyn_5scores[[name]], 0.95))
  syn_rows <- rownames(combosyn_5scores[combosyn_5scores[[name]] >= quantile(combosyn_5scores[[name]], 0.95),])
  an_rows <-rownames(combosyn_5scores[combosyn_5scores[[name]] < quantile(combosyn_5scores[[name]], 2/3),])
  
  new_col <- paste(name,"_label")
  combosyn_5scores <- combosyn_5scores  %>% 
    mutate(!!new_col := case_when(
      rownames(combosyn_5scores) %in%  syn_rows ~ "synergistic",
      rownames(combosyn_5scores) %in%  an_rows ~ "antagonistic",
    ))
}


combosyn_5scores <- assigning(name="HSA")
combosyn_5scores <- assigning(name="Bliss")
combosyn_5scores <- assigning(name="Loewe")
combosyn_5scores <- assigning(name="ZIP")
combosyn_5scores <- assigning(name="COMBO_SCORE")

upframe <- combosyn_5scores%>% 
  dplyr::select("HSA _label","Bliss _label","Loewe _label",
                "ZIP _label","COMBO_SCORE _label")
names(upframe) = c("HSA","Bliss","Loewe","ZIP","COMBO_SCORE")

#upsetplot for synergistic pair
temp_syn_up <- upframe %>% mutate(across(.fns = ~replace(., . !=  "synergistic" , 0)))
temp_syn_up <- temp_syn_up %>% mutate(across(.fns = ~replace(., . ==  "synergistic" , 1)))
temp_syn_up[is.na(temp_syn_up)] <- 0
temp_syn_up <- temp_syn_up %>%
  mutate_if(is.character, as.numeric)
fig1d <- UpSetR::upset(temp_syn_up, sets=c("COMBO_SCORE","HSA","Bliss","ZIP","Loewe"),
              order.by = c("freq","degree"), decreasing = c(TRUE,TRUE),
              main.bar.color = "coral2",mainbar.y.label="Overlap of Synergistic Pairs",
              keep.order = TRUE)
fig1d

temp_syn_an <- upframe %>% mutate(across(.fns = ~replace(., . !=  "antagonistic" , 0)))
temp_syn_an <- temp_syn_an %>% mutate(across(.fns = ~replace(., . ==  "antagonistic" , 1)))
temp_syn_an[is.na(temp_syn_an)] <- 0
temp_syn_an <- temp_syn_an %>%
  mutate_if(is.character, as.numeric)
fig1d_2 <- UpSetR::upset(temp_syn_an, sets=c("COMBO_SCORE","HSA","Bliss","ZIP","Loewe"),
                       order.by = c("freq","degree"), decreasing = c(TRUE,TRUE),
                       main.bar.color = "aquamarine3",mainbar.y.label="Overlap of Antagonistic Pairs",
                       keep.order = TRUE)
fig1d_2

# fig1b
names(combosyn_5scores) = c(names(combosyn_5scores)[1:13], c("HSA_label","Bliss_label","Loewe_label",
                                                            "ZIP_label","COMBO_SCORE_label"))
combosyn_5scores <- combosyn_5scores %>% group_by(pubchemID_x,pubchemID_y) %>% mutate(pairs_id = cur_group_id())

jaccard_mtx <- function(name) {
  HSA_wide <- combosyn_5scores %>% 
    filter(Bliss_label == "synergistic") %>% 
    dcast(., Panel.Name~pairs_id)  %>%
    column_to_rownames("Panel.Name") %>% 
    mutate_if(is.numeric, ~1 * (. > 0)) %>%
    mutate_if(is.character, as.numeric)
  hsa.dj <- vegdist(HSA_wide, "jac", binary=TRUE)
  hsa.dj
}
hsa.dj <- jaccard_mtx(HSA_label)
bliss.dj <- jaccard_mtx(Bliss_label)
loewe.dj <- jaccard_mtx(Loewe_label)
zip.dj <- jaccard_mtx(ZIP_label)
combo.dj <- jaccard_mtx(COMBO_SCORE_label)
avg.dj.ls<- list(hsa.dj,bliss.dj,loewe.dj,zip.dj,combo.dj)
avg.dj <- reduce(avg.dj.ls,`+`) / length(avg.dj.ls)
# average metrics
M <- 1-as.matrix(avg.dj)
M <- M[, c("Prostate","Breast","Leukemia","CNS","Colon","Ovarian","Renal","Melanoma","Non-Small Cell Lung")]
M <- M[c("Prostate","Breast","Leukemia","CNS","Colon","Ovarian","Renal","Melanoma","Non-Small Cell Lung"),]
M[lower.tri(M)] <- NA
diag(M) <- NA

pheatmap(M, cluster_rows=F, cluster_cols=F, na_col="white",
         border_color ="white",color=scales::seq_gradient_pal("white", "tomato3", "Lab")(seq(0,1,length.out=30)))

combosyn_5scores %>% 
  ggplot(aes(x=fct_infreq(Panel.Name))) + geom_bar(stat="count",fill = "#FF6666") +
  coord_flip() + labs(x= "Tissue", y= "Set Size") +
  theme_classic()

#fig 2a -  python to calculate similarity
#write.table(combosyn_5scores, file = "combosyn_5scores.csv",sep=",")
combosyn_5scores <- read.csv(paste0(data_dir, "/combosyn_5scores_dbID.csv"))

combosyn_5scores <- assigning(name="HSA")
combosyn_5scores <- assigning(name="Bliss")
combosyn_5scores <- assigning(name="Loewe")
combosyn_5scores <- assigning(name="ZIP")
combosyn_5scores <- assigning(name="COMBO_SCORE")
#
chem_melt <- combosyn_5scores_pths %>% 
  dplyr::select("HSA_label","Bliss_label","Loewe_label",
                "ZIP_label","COMBO_SCORE_label","chemical_similarity") %>% 
  rename(., "HSA"="HSA_label","Bliss"="Bliss_label","Loewe"="Loewe_label","ZIP"="ZIP_label","ALMANAC SCORE"="COMBO_SCORE_label") %>% 
  melt(., id.vars="chemical_similarity", variable.name="metrics",value.name = "type")
chem_melt <- chem_melt %>% mutate(type = case_when(
  type %in% c("") ~ "Neither",
  type %in% c("synergistic") ~ "Synergistic",
  type %in% c("antagonistic") ~ "Antagonistic"))

ggplot(chem_melt, aes(chemical_similarity, colour = type,linetype=metrics)) +
  stat_ecdf(geom = "step") + scale_color_manual(values = c("#91D1C2B2",'grey',"#E64B35B2")) +
  theme_classic()


# get pathway for each drug
# kegg
drug_bag <- unique(c(combosyn_5scores$compound1_x,combosyn_5scores$compound1_y))
dpi <- read.csv("../synergyy/results/proessed_dpi_db.csv", header=FALSE, sep=",")
names(dpi) <- dpi[1,]
dpi <- dpi[-1,]
rownames(dpi) <- dpi[,1]
dpi[,1] <- NULL

#KEGG 
kegg_pths = list()
for(mol in drug_bag){
  print(as.character(mol))
  v <- as.character(mol)
  #only include drug that have target genes > 2
  if (length(rownames(dpi[dpi[[v]] ==1,])) >= 2){
    # kk <- enrichKEGG(gene   = rownames(dpi[dpi[[v]] ==1,]), #
    #                  organism     = 'hsa',
    #                  minGSSize = 5,
    #                  pvalueCutoff = 0.05)
    #print(kk$GeneRatio/kk$BgRatio)
    kk <-  enrichPathway(gene   = rownames(dpi[dpi[[v]] ==1,]))
    kegg_pths[[v]] <- kk$Description
  }
  
}

DICE <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  J = intersection/union
  S = 2*J/(1+J)
  return (S)
}

combosyn_5scores_pths <-combosyn_5scores %>% filter(compound1_x %in% names(kegg_pths) &
                              compound1_y %in% names(kegg_pths))

combosyn_5scores_pths$pth_similarity <- NA
for(row in 1:nrow(combosyn_5scores_pths)) {
    d1 <- as.character(combosyn_5scores_pths[row, "compound1_x"])
    d2  <- as.character(combosyn_5scores_pths[row, "compound1_y"])
    k1 <- kegg_pths[[d1]]
    k2 <-  kegg_pths[[d2]]
    combosyn_5scores_pths[row, "pth_similarity"] <- DICE(k1,k2)
}

pths_melt <- combosyn_5scores_pths %>% 
  dplyr::select("HSA_label","Bliss_label","Loewe_label",
                "ZIP_label","COMBO_SCORE_label","pth_similarity") %>% 
  rename(., "HSA"="HSA_label","Bliss"="Bliss_label","Loewe"="Loewe_label","ZIP"="ZIP_label","ALMANAC SCORE"="COMBO_SCORE_label") %>% 
  melt(., id.vars="pth_similarity", variable.name="metrics",value.name = "type")
pths_melt <- pths_melt %>% mutate(type = case_when(
  type %in% c("") ~ "Neither",
  type %in% c("synergistic") ~ "Synergistic",
  type %in% c("antagonistic") ~ "Antagonistic"))

ggplot(pths_melt, aes(pth_similarity, colour = type,linetype=metrics)) +
  stat_ecdf(geom = "step") + scale_color_manual(values = c("#91D1C2B2",'grey',"#E64B35B2")) +
  theme_classic()












