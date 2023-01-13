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
library(OmnipathR)
library(igraph)
library(ggraph)
library(janitor)

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
  dplyr::select("HSA _label","Bliss _label","Loewe _label",
                "ZIP _label","COMBO_SCORE _label","chemical_similarity") %>% 
  #rename(., "HSA"="HSA _label","Bliss"="Bliss _label","Loewe"="Loewe _label","ZIP"="ZIP _label","ALMANAC SCORE"="COMBO_SCORE _label") %>% 
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
    kk <- enrichKEGG(gene   = rownames(dpi[dpi[[v]] ==1,]), #
                     organism     = 'hsa',
                     minGSSize = 5,
                     pvalueCutoff = 0.05)
    #print(kk$GeneRatio/kk$BgRatio)
    #kk <-  enrichPathway(gene   = rownames(dpi[dpi[[v]] ==1,]))
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
  dplyr::select("HSA _label","Bliss _label","Loewe _label",
                "ZIP _label","COMBO_SCORE _label","pth_similarity") %>% 
  #rename(., "HSA"="HSA _label") %>% 
  melt(., id.vars="pth_similarity", variable.name="metrics",value.name = "type")
pths_melt <- pths_melt %>% mutate(type = case_when(
  type %in% c("") ~ "Neither",
  type %in% c("synergistic") ~ "Synergistic",
  type %in% c("antagonistic") ~ "Antagonistic"))

ggplot(pths_melt, aes(pth_similarity, colour = type,linetype=metrics)) +
  stat_ecdf(geom = "step") + scale_color_manual(values = c("#91D1C2B2",'grey',"#E64B35B2")) +
  theme_classic()

#flatten list
melt_kegg_pths <- kegg_pths %>% do.call(rbind, .) %>% melt(.) %>% 
  dplyr::select("Var1","value",d=Var1,pth=value) %>% distinct()

# kegg,
r <- combosyn_5scores_pths %>% rename("HSA_label"="HSA _label","Bliss_label"="Bliss _label","Loewe_label"="Loewe _label","ZIP_label"="ZIP _label","COMBO_SCORE_label"="COMBO_SCORE _label")
#
subset <- r %>% dplyr::select("compound1_x","compound1_y","HSA_label")  %>% 
  group_by(compound1_x,compound1_y) %>% 
  count(HSA_label, sort = TRUE) %>% dcast(.,compound1_x+compound1_y~HSA_label)  %>% 
  dplyr::select("compound1_x","compound1_y","antagonistic","synergistic") 
subset[is.na(subset)] <- 0

subset_pth <- subset %>% 
  left_join(., melt_kegg_pths,by = c("compound1_x" = "d")) %>% 
  left_join(., melt_kegg_pths,by = c("compound1_y" = "d")) %>% 
  drop_na()
  #remove rows with na

#This is background:sum(subset$antagonistic) = 133228 records; sum(subset$synergistiv) = 11571 records
subset_pth_sum <- subset_pth %>% group_by(pth.x,pth.y) %>% summarise(synergistic_sum=sum(synergistic),
                                                                     antagonistic_sum=sum(antagonistic))
#Toy sample,
# 144799 records of which 11571 are synergistic,
# a pathway(combo) with 1149+3359=4508 records,
# including 1149 of the 11571 synergy records
##                 syn             anta
## pathwayYes    1149           3359
## pathwayNo    11571-1149                             144799-4508-(11571-1149) 

# dat <- matrix(c(1149, 10422, 3359, 129869), nrow = 2, ncol = 2)
# rownames(dat) <- c("Mutated gene", "No mutated gene")
# colnames(dat) <- c("Cancer", "Normal")
# fisher.test(dat)

myFUN<- function(x) {
  py_syn = as.numeric(x["synergistic_sum"])
  py_anta = as.numeric(x["antagonistic_sum"])
  dat <- matrix(c(py_syn, (11571-py_syn), py_anta, (144799-(py_syn+py_anta)-(11571-py_syn))), nrow = 2, ncol = 2)
  t <- fisher.test(dat)
  # x["low"] <- t$conf.int[1]
  # x["high"] <- t$conf.int[2]
  # x["est"] <- t$estimate
  c(t$conf.int[1],t$conf.int[2],t$estimate)
  }
or_output <- as.data.frame(t(apply(subset_pth_sum,1,myFUN)))
subset_pth_sum$low_HSA <- or_output$V1
subset_pth_sum$high_HSA <- or_output$V2
subset_pth_sum$odds_HSA <- or_output$"odds ratio"

# write.table(subset_pth_sum, file = "or_HSA.csv",sep=",",row.names=FALSE)

or_hsa <- read.csv(paste0(data_dir, "/or_HSA.csv"))
or_bliss <- read.csv(paste0(data_dir, "/or_Bliss.csv"))
or_loewe <- read.csv(paste0(data_dir, "/or_Loewe.csv"))
or_zip <- read.csv(paste0(data_dir, "/or_ZIP.csv"))
or_combo <- read.csv(paste0(data_dir, "/or_COMBO_SCORE.csv"))

MergedDF <- Reduce(function(x,y) merge(x,y,by= c("pth.x","pth.y"),all=TRUE),
                   list(or_hsa,or_bliss,or_loewe,or_zip,or_combo))
#at least one synergy metric with an Odd’s Ratio lower confidence interval above 1.5 
#and an Odd’s Ratio 373higher confidence interval lower than 1. 


MergedDF_t <- MergedDF %>% dplyr::select(starts_with("pth") | starts_with("low") | starts_with("high") | starts_with("odds")) %>% 
  filter_at(vars(starts_with("low")), any_vars(. >= 1.5)) %>%  
  filter_at(vars(starts_with("high")), any_vars(. < 1)) %>%
  dplyr::select(starts_with("pth") | starts_with("odds")) %>% 
  pivot_longer(-c("pth.x","pth.y"), names_to = "group", values_to = "y")

ggplot(MergedDF_t, aes(y, colour =group)) +
  stat_ecdf(geom = "step")+
  xlim(0, 5) + scale_x_continuous(limits=c(0,5), breaks=seq(0,5, by = 1))+
  labs(x= "Odds ratio", y= "Frequency") +
  theme_classic() 

MergedDF_s <- MergedDF %>% dplyr::select(starts_with("pth") | starts_with("low") | starts_with("high") | starts_with("odds")) %>% 
  filter_at(vars(starts_with("low")), any_vars(. >= 1.5)) %>%  
  filter_at(vars(starts_with("high")), any_vars(. < 1))
  
MergedDF_s <- MergedDF_across %>% filter(pth.x=="Rap1 signaling pathway",pth.y=="Phagosome") 
df <- data.frame(yAxis = c("ZIP","Loewe","HSA","ALMANAC","Bliss"), 
                 boxOdds = c(MergedDF_s$odds_ZIP,MergedDF_s$odds_Loewe,MergedDF_s$odds_HSA,MergedDF_s$odds_COMBO_SCORE,MergedDF_s$odds_Bliss), 
                 boxCILow = c(MergedDF_s$low_ZIP,MergedDF_s$low_Loewe,MergedDF_s$low_HSA,MergedDF_s$low_COMBO_SCORE,MergedDF_s$low_Bliss), 
                 boxCIHigh = c(MergedDF_s$high_ZIP,MergedDF_s$high_Loewe,MergedDF_s$high_HSA,MergedDF_s$high_COMBO_SCORE,MergedDF_s$high_Bliss)
)

fig2c <- ggplot(df, aes(x = boxOdds, y = yAxis)) +
            geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
            geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height =
                             .2, color = "gray50") +
            geom_point(size = 3.5, color = "orange")+
            labs(x= "Odds ratio", y= "") + xlim(0, 13) +
              theme_classic()
fig2c

## Do we also find pathway combination that are irrelvant with metrics calculation
MergedDF_across <- MergedDF %>% dplyr::select(starts_with("pth") | starts_with("low") | starts_with("high") | starts_with("odds")) %>% 
  filter_at(vars(starts_with("low")), all_vars(. >= 1.5))

#Figure2D, https://workflows.omnipathdb.org/drug-targets.html#initialise-omnipath-database
#https://r.omnipathdb.org/reference/OmnipathR.html
combosyn_5scores <- read.csv(paste0(data_dir, "/combosyn_5scores_dbID.csv"))
dpi <- read.csv("processed_dpi.csv", header=FALSE, sep=",")
names(dpi) <- dpi[1,]
dpi <- dpi[-1,]
rownames(dpi) <- dpi[,1]
dpi[,1] <- NULL

# Download protein-protein interactions
interactions <- import_omnipath_interactions()
OPI_g <- interaction_graph(interactions = interactions )

#ncbi2symbol
ncbi2symbol <- read.csv("ncbi2symbol.txt", header=TRUE, sep="\t")
colnames(ncbi2symbol) <- c("symbol","id")
ncbi2symbol$id <- as.character(ncbi2symbol$id)
#same target where drug combs had at least one drug target in common

# same pathway
# shortest_paths(OPI_g, from = "STAT3", to = "CTTN",output = 'epath')$epath[[1]]
# shortest_paths(OPI_g, from = "CTTN", to = "STAT3",output = 'epath')$epath[[1]]

#parallel pathway: drug targets were upregulated or downregulated by the same gene
# as.character(incident(OPI_g, "STAT3", mode=c("in")))
# as.character(incident(OPI_g, "CTTN", mode=c("in"))) 

#
combosyn_5scores_net <- combosyn_5scores %>% distinct(compound1_x, compound1_y, .keep_all = FALSE)
#移除没有target的药物
drug_bag_in_dpi <- unique(colnames(dpi))
combosyn_5scores_net <-combosyn_5scores_net %>% filter(compound1_x %in% drug_bag_in_dpi &
                                                      compound1_y %in% drug_bag_in_dpi)

combosyn_5scores_net$network_class <- NA
for(row in 1:nrow(combosyn_5scores_net)) {
  d1 <- as.character(combosyn_5scores_net[row, "compound1_x"])
  d2  <- as.character(combosyn_5scores_net[row, "compound1_y"])
  
  tg1 <- dpi %>% dplyr::select(d1)  %>% `colnames<-`(c("drug_name")) %>% 
    mutate(index = rownames(.)) %>% 
    filter(., drug_name == 1) %>% dplyr::select("index")
  tg2 <- dpi %>% dplyr::select(d2)  %>% `colnames<-`(c("drug_name")) %>% 
    mutate(index = rownames(.)) %>% 
    filter(., drug_name == 1) %>% dplyr::select("index")
  
  #首先判断是否有相同的target
  if(nrow(intersect(tg1,tg2))>0){
    combosyn_5scores_net[row, "network_class"] <- "same target"
  }else{
    #判断是否有属于同一通路
   
    # 遍历所有的基因
    empty_list <- c()
    graph_vertex <- V(OPI_g)$name
    for(i in 1:nrow(tg1)){
      for(j in 1:nrow(tg2)){
        # 轉換成gene 名
        g1 <- as.character(ncbi2symbol %>% dplyr::filter(id %in% tg1[i,]) %>% dplyr::select(symbol))
        g2 <- as.character(ncbi2symbol %>% dplyr::filter(id %in% tg2[j,]) %>% dplyr::select(symbol))
        
        if (g1 %in% graph_vertex && g2 %in% graph_vertex){
          sp1 <- shortest_paths(OPI_g, from = g1, to = g2,output = 'epath')$epath[[1]]
          sp2 <- shortest_paths(OPI_g, from = g2, to = g1,output = 'epath')$epath[[1]]
          
          if(length(sp1)>0 || length(sp2)>0){
            empty_list <- rbind(empty_list,"same pathway")
          }else{
            # 判断是否属于平行通路
            in_g1 <- as.character(neighbors(OPI_g, g1, mode=c("in")))
            in_g2 <- as.character(neighbors(OPI_g, g2, mode=c("in")))
            in_ <- length(intersect(in_g1,in_g2))
            if(!is.null(in_)){
              empty_list <- rbind(empty_list,"parallel pathway")
            }else{
              empty_list <- rbind(empty_list,"unidenfitied")
            }
          }
          
        }else{# 在graph里找不到
          empty_list <- rbind(empty_list,"notingraph")
        }


      
      }
    }
    combosyn_5scores_net[row, "network_class"] <- empty_list[1]
  }
}

write.table(combosyn_5scores_net, file = "combosyn_network_class.csv",sep=",")
combosyn_network_class <- read.csv(paste0(data_dir, "/combosyn_network_class.csv"))

combosyn_5scores_network <- left_join(combosyn_5scores, combosyn_network_class,
                                    by = c("compound1_x" = "compound1_x", "compound1_y" = "compound1_y"))
combosyn_5scores_network <- combosyn_5scores_network %>%
  mutate(network_class = if_else(is.na(network_class), "notingraph",network_class))
#Here we want to test whether syn/non differ between same target and others for each synergy metric
# Start from ALMNAC score, 
# same target/same pathway/parallel pathway v.s. unidenfitied 
# combosyn_5scores_network <- combosyn_5scores_network %>% filter(network_class %in%  c("unidentified","same pathway",
#                                                         "same target","parallel pathway"))



# #试比较same target组HSA是否比其他组hsa高？
# 

#create an empty df to store w test value
w_df <-c("yAxis","boxOdds","boxCILow","boxCIHigh")

w_test <- function(score, score_label, target_class){

  a <- combosyn_5scores_network %>% 
    mutate(network_class_wx = case_when(
      network_class == "notingraph" ~ "others",
      # network_class == "unidenfitied " ~ "others",
      # network_class != target_class ~ "others",
      network_class == target_class ~ target_class))
  
  a_ <- a %>% filter(!! as.name(score_label) %in%  c("synergistic","antagonistic"))
  # a_ <- a
  x <- a_ %>% filter(network_class_wx %in%  c(!!target_class)) %>% select(all_of(score))
  y <- a_ %>% filter(network_class_wx %in%  c("others")) %>% select(all_of(score))
  Nj  <- c(nrow(x), nrow(y))
  wIndDf <- rbind(x,y)
  wIndDf$IV <- factor(rep(1:2, Nj), labels=LETTERS[1:2])
  colnames(wIndDf) <- c("v","IV")
  w <- wilcox.test(v ~ IV, conf.int=TRUE, data=wIndDf,alternative="two.sided",tol.root = 1e-4)
  w_df <-rbind(w_df,c(score,w$estimate,w$conf.int[1],w$conf.int[2]))
  # w_df$yAxis <- "HSA"
  # w_df$boxOdds <- w$estimate
  # w_df$boxCILow <- w$conf.int[1]
  # w_df$boxCIHigh <- w$conf.int[2]
  
}
w_df <- w_test("COMBO_SCORE","COMBO_SCORE _label","same target")
w_df <- w_test("ZIP","ZIP _label","same target")
w_df <- w_test("Loewe","Loewe _label","same target")
w_df <- w_test("HSA","HSA _label","same target")
w_df <- w_test("Bliss","Bliss _label","same target")
#       
w_df <- data.frame(w_df) %>% row_to_names(row_number = 1)
# 

w_df_ <- w_df %>%
  mutate_if(is.character, as.numeric)
w_df_$yAxis <- c("COMBO_SCORE","ZIP","Loewe","HSA","Bliss")

fig2d <- ggplot(w_df_, aes(x = boxOdds, y = yAxis)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height =
                   .2, color = "gray50") +
  geom_point(size = 2, color = "orange")+
  labs(x= "Wilcox Location Shift", y= "") + xlim(-1, 6) +
  theme_classic()
fig2d
  
  
  
  
  
  
  


