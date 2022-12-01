library(synergyfinder)
library(data.table)
library(dplyr)
# library(classyfireR)
# library(ggplot2)
# library(forcats)

data_dir = "/Users/chengqi_xu/Documents/Elemento lab/clinicaltrials.gov"
# Load the ALMANAC data.
almanac_dat = read.csv(paste0(data_dir, "/ComboDrugGrowth_Nov2017.csv"))
# Convert it to a data.table. 
setDT(almanac_dat)

# Calculate the CORRECT PercentGrowth. T_T
almanac_dat$PERCENTGROWTHCORRECTED = 100 * (almanac_dat$TESTVALUE - almanac_dat$TZVALUE) / (almanac_dat$CONTROLVALUE - almanac_dat$TZVALUE)
head(almanac_dat)

mylookup_tibble = read.csv(paste0(data_dir, "/ComboCompoundNames_small.txt"),sep='\t',)
joined_tibble <- left_join(almanac_dat, mylookup_tibble, 
                           by = c("NSC1" = "IX"))
joined_tibble2 <- left_join(joined_tibble, mylookup_tibble, 
                           by = c("NSC2" = "IX"))
# create a block_id for groupby
joined_tibble2 <- joined_tibble2 %>% group_by(Drug.x,Drug.y,CELLNAME) %>% mutate(block_id = cur_group_id())

screening_data <- joined_tibble2 %>% 
  dplyr::select(block_id,Drug.x,Drug.y,CONC1, CONC2, PERCENTGROWTHCORRECTED, CONCUNIT1,CONCUNIT2) %>% 
  rename( "drug_row" = "Drug.x", 
          "drug_col" = "Drug.y",
          "conc_r" = "CONC1",
          "conc_c" = "CONC2", 
          "response"="PERCENTGROWTHCORRECTED", 
          "conc_r_unit"="CONCUNIT1",
          "conc_c_unit"="CONCUNIT2")

screening_data$conc_r <- screening_data$conc_r * 10^9
screening_data$conc_c <- screening_data$conc_c * 10^9
screening_data$conc_r_unit <- "nM"
screening_data$conc_c_unit <- "nM"

# 60 cancer cell lines
# Fig1a. drugbank -- superclass
batch_classifire = read.csv(paste0(data_dir, "/classyfire_.csv"))
setDT(batch_classifire)
g1 <- batch_classifire %>% mutate(Superclass = Superclass %>% fct_infreq()  %>% fct_rev()) %>% 
  ggplot(aes(x=Superclass,fill=Superclass)) +
  geom_bar(stat="count") + coord_flip() +
  scale_fill_manual(values=scales::seq_gradient_pal("white", "dodgerblue4", "Lab")(seq(0,1,length.out=15))) + 
  labs(x= "Drug Class", y= "Set Size") +
  theme_classic()



data("mathews_screening_data")
head(mathews_screening_data)

#fig1b
#需要增加monotherapy的数据
screening_data_mono <- screening_data[is.na(screening_data$drug_col),]

test <- screening_data %>% filter(block_id == 282352)


res <- synergyfinder::ReshapeData(
  data = test,
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

#fig1c

