

####Microbiome analysis#####


################### crypotpoc study ###################

#install.packages("tidyverse")

library(tidyverse)
library(ggpubr)
library(ggsci)
library(dplyr)
library(rstatix)
library(readxl)
library(ampvis2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(phyloseq)
library(vegan)
library(DESeq2)
library(EnhancedVolcano)
library(ggrepel)
library(metagMisc)
library(RColorBrewer)
library(ComplexHeatmap)
library(metagenomeSeq)
library(microbiomeMarker)
library(microbiome)
library(export)
library(rabuplot)

#Coulour pallettes
col_fil <- pal_jco("default")(10)

col_scale <- scale_color_jco()

#Set working directory to script directory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Load phyloseq objects


source("Phyloseq_import_16S.R")

source("Phyloseq_CSS_import_16S.R")


## recategorizing variables##

sample_data(PSB)$diarrhoea_category<- revalue(sample_data(PSB)$diarrhoea,c("AD"="AD", "ProD"="ProPD", "PD"="ProPD","Cont"="Control"))  ## categorizing diarrhea  into two categories 



PSB.CSS
sample_data(PSB.CSS)$diarrhoea_category<- revalue(sample_data(PSB.CSS)$diarrhoea,c("AD"="AD", "ProD"="ProPD", "PD"="ProPD","Cont"="Control"))  # categorizing diarrhea  into two catagories 

### subset Case_control  samples only by removing the follow up samples ######## 

PSOBT<- subset_samples(PSB,cc_substudy!="0") ## this are only case control samples follow up samples removed. 
PSOBT

PSOBT.CSS <-subset_samples(PSB.CSS,cc_substudy!="0")  ## follow up samples removed 
PSOBT.CSS


PSOB<- subset_samples(PSOBT,control_fct!="NA") ## this are only case control samples follow up samples removed. 
PSOB

sample_data(PSOB)$month<-as.factor(sample_data(PSOB)$month)## to change the numric variable to factor so that it will be easy to convert categorical variable blow

sample_data(PSOB)$season<-revalue(sample_data(PSOB)$month,c("2"="dry_season","3"="dry_season","4"="dry_season","5"="wet_season","6"="wet_season", "7"="wet_season","8"="wet_season", "9"="wet_season","10"="wet_season","11"="dry_season", "12"="dry_season", "13"="dry_season", "14"="dry_season","15"="dry_season","16"="dry_season","17"="wet_season","18"="wet_season","19"="wet_season")) ## classifying the seasons to category 


PSOB.CSS <-subset_samples(PSOBT.CSS,control_fct!="NA")  ## follow up samples removed 
PSOB.CSS


sample_data(PSOB.CSS)$month<-as.factor(sample_data(PSOB.CSS)$month)## to change the numeric variable to factor 

sample_data(PSOB.CSS)$season<-revalue(sample_data(PSOB.CSS)$month,c("2"="dry_season","3"="dry_season","4"="dry_season","5"="wet_season","6"="wet_season", "7"="wet_season","8"="wet_season", "9"="wet_season","10"="wet_season","11"="dry_season", "12"="dry_season", "13"="dry_season", "14"="dry_season","15"="dry_season","16"="dry_season","17"="wet_season","18"="wet_season","19"="wet_season")) ## classifying the seaseons



PSBJ <- prune_samples(sample_sums(PSOB) >=5000, PSOB) #### remove samples below 5000 read per samples

PSBJ


PSBJ.CSS <- prune_samples(sample_sums(PSOB) >=5000, PSOB.CSS) #### remove samples below 5000 read per samples  in the CSS

PSBJ.CSS


#### remove zOTUs   found  in only 1% of the samples ##########


PSBJ.filtered <- phyloseq_filter_prevalence(PSBJ, prev.trh=0.01) ## remove zOTU found only 1% of samples



PSBJ.filtered #### this used phyloseq for the further all  analysis alhadiversity
 
PSBJ.CSS.filtered <- phyloseq_filter_prevalence(PSBJ.CSS, prev.trh = 0.01) ## remove zOTU found only 1% of samples

PSBJ.CSS.filtered   #### This the  phyloseq object used for  betaadversity analysis



############################## Diarrhea microbiome study ################

#######################Basic stats####################################



#######################Bar charts#######################


###genus relative  abundance ##########


##age_stratum relative abundance

ps0 <- merge_samples(ps.genus.rel, "age_stratum")

ps0 <- transform_sample_counts(ps0, function(x) x / sum(x)) #collapsed otu table to relative abundance

#Create melted dataframe
df <- psmelt(ps0)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)]," unknown")

#Arrange samples by mean abundance
top <- df %>%
  group_by(Sample, tax) %>%
  dplyr::summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
#Show top
top

top10 <- top$tax[1:45] #Select  most abundant genera
df0 <- df %>% 
  mutate(tax = fct_other(tax, c(as.matrix(top10)))) #Combine all other entries as "other"

df0 <- df0[order(df0$Sample, decreasing = TRUE), ]

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = c("0-5m", "6-11m","12-23m","24-60m"))

barplot.age_stratum_genus<- ggplot(df0, aes(Sample, Abundance, fill = fct_reorder(tax, Abundance, .fun = mean))) + 
  geom_col(width = 0.8) +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_fil,5),name="Genus")+
  theme_classic() +
  theme(text = element_text(size = 10, colour = "Black"),
        axis.line=element_line(size=0.5),
        #panel.border = element_blank(),
        axis.text=element_text(size = 10, colour = "Black"),
        axis.ticks=element_line(size=1, colour = "Black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(angle = 0, size=12, face = "bold"),
        #legend.position = "none"
  ) +
  ylab("Relative abundance")
#ggtitle("age_stratum")#


barplot.age_stratum_genus


# #####Top 20 taxa by age groups stratified by cases and controls#### 

PSBJ.filtered@sam_data$age_stratum <- factor(PSBJ.filtered@sam_data$age_stratum, levels=c("0-5m", "6-11m","12-23m", "24-60m"))


library(rabuplot)### to plot the top 20 taxa

Age_groupestop20taxa <- rabuplot(PSBJ.filtered,
                                 violin = TRUE, predictor = "age_stratum",N_taxa=20,By_median =FALSE, type="genus",no_other_type=TRUE, facet_wrap= "control_fct",colors=col_fil[1:4]) 

Age_groupestop20taxa



################# Diarrhea category relative abundance###################

ps.genus <- tax_glom(PSBJ.filtered, "Genus", NArm = FALSE) #Collapse to genus level the Phyloseq object 
ps.genus.rel <- transform_sample_counts(ps.genus, function(x) x / sum(x)) #collapsed otu table to relative abundance

ps0 <- merge_samples(ps.genus.rel,"diarrhoea_category")

ps0 <- transform_sample_counts(ps0, function(x) x / sum(x))

#Create melted dataframe
df <- psmelt(ps0)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)]," unknown")

#Arrange samples by mean abundance
top <- df %>%
  group_by(Sample, tax) %>%
  dplyr::summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
#Show top
top

top10 <- top$tax[1:40] #Select the top  most abundant genera
df0 <- df %>% 
  mutate(tax = fct_other(tax, c(as.matrix(top10)))) #Combine all other entries as "other"

df0 <- df0[order(df0$Sample, decreasing = TRUE), ]


df0$Sample <- factor(df0$Sample, levels = c("Control","AD", "ProPD"), labels= c("Control","AD", "ProPD"))


barplot.diarrhoea_duration_genus <- ggplot(df0[!is.na(df0$Sample),],aes(Sample, Abundance, fill = fct_reorder(tax, Abundance, .fun = mean))) + 
  geom_col(width = 0.8) +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_fil,5),name="Genus")+
  theme_classic() +
  theme(text = element_text(size = 10, colour = "Black"),
        axis.line=element_line(size=0.5),
        #panel.border = element_blank(),
        axis.text=element_text(size = 10, colour = "Black"),
        axis.ticks=element_line(size=1, colour = "Black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(angle = 0, size=12, face = "bold"),
        #legend.position = "none"
  ) +
  ylab("Relative abundance") 
#ggtitle("diarrhoea_duration")#


barplot.diarrhoea_duration_genus



#############################Alpha diversity##########################


#Set desired alpha diversity metrics for all alphadiversity measures## 


alha_met <- c("Observed", "Shannon")# ,"Chao1","ACE", "InvSimpson", "Fisher")


#Set sign argumetns
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


##########Cases_control alpha diversity############

my_comparisons = list( c("Case", "Control")) #All comparisons



bac.alpha.Case_Control <- plot_richness(PSBJ.filtered, measures=alha_met , x="control_fct", color="control_fct") + 
  geom_boxplot(alpha=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = FALSE, symnum.args = symnum.args) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank()
  ) +
  geom_point(size = 2) +
  scale_color_manual(values=col_fil) +
  labs(y = "Alpha diversity") 
#ggtitle("Alpha diversity_Case_Control")##

bac.alpha.Case_Control


####Duration Diarrhoea alpha measures  

PSBJ.filtered@sam_data$diarrhoea_category<- factor(PSBJ.filtered@sam_data$diarrhoea_category, levels=c("Control", "AD", "ProPD"))


my_comparisons = list(c("Control","AD"), c("Control", "ProPD"), c("AD", "ProPD")) #All comparisons


PSB.sub_alpha_duration <- subset_samples(PSBJ.filtered, diarrhoea_category != "NA")##### to remove the NA in  "diarrhoea" variable####

bac.alpha.duration_diarrhoea <- plot_richness(PSB.sub_alpha_duration, measures=alha_met , x="diarrhoea_category", color="diarrhoea_category") + 
  geom_boxplot(alpha=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = FALSE, symnum.args = symnum.args) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank()
  ) +
  geom_point(size = 2) +
  scale_color_manual(values=c("Control"="#EFC000","AD"="#0173C2","ProPD"= "#868686")) +
  labs(y = "Alpha diversity") 
#ggtitle("Alpha diversity_duration of diarrhe")##

bac.alpha.duration_diarrhoea



######Age_stratum alpha diversity measure ######

PSBJ.filtered@sam_data$age_stratum <- factor(PSBJ.filtered@sam_data$age_stratum, levels=c("0-5m", "6-11m","12-23m", "24-60m"))

my_comparisons = list( c("0-5m", "6-11m"),c("0-5m","12-23m"),c("0-5m","24-60m"),c("6-11m","12-23m"), c("6-11m","24-60m"),c("12-23m","24-60m")) #All comparisons

 
alpha_age_stratum<- plot_richness(PSBJ.filtered, measures=alha_met , x="age_stratum", color="age_stratum") + 
  geom_boxplot(alpha=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = FALSE, symnum.args = symnum.args) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank()
  ) +
  geom_point(size = 2) +
  scale_color_manual(values=col_fil) +
  labs(y = "Alpha diversity") 
#ggtitle("Alpha diversity_age_stratum)##
alpha_age_stratum



##### Alpha diversity cases and controls stratified by age stratum######

####subset controls#####

PSB.sub_alpha_controls_only<- subset_samples(PSBJ.filtered,control_fct!= "Case") ##  subset controls only  


PSB.sub_alpha_controls_only@sam_data$age_stratum <- factor(PSB.sub_alpha_controls_only@sam_data$age_stratum, levels=c("0-5m", "6-11m","12-23m", "24-60m"))

my_comparisons = list( c("0-5m", "6-11m"),c("0-5m","12-23m"),c("0-5m","24-60m"),c("6-11m","12-23m"), c("6-11m","24-60m"),c("12-23m","24-60m")) #All comparisons

 
alpha_age_stratum_control<- plot_richness(PSB.sub_alpha_controls_only, measures=alha_met , x="age_stratum", color="age_stratum") + 
  geom_boxplot(alpha=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = FALSE, symnum.args = symnum.args) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank()
  ) +
  geom_point(size = 2) +
  scale_color_manual(values=col_fil) +
  labs(y = "Alpha diversity") 
#ggtitle("Alpha diversity_age_stratum_cotrol")##
alpha_age_stratum_control



### ##############Subset Cases###############################

PSB.sub_alpha_case_only<- subset_samples(PSBJ.filtered,control_fct!= "Control") ##subset cases


PSB.sub_alpha_case_only@sam_data$age_stratum <- factor(PSB.sub_alpha_case_only@sam_data$age_stratum, levels=c("0-5m", "6-11m","12-23m", "24-60m"))

my_comparisons = list( c("0-5m", "6-11m"),c("0-5m","12-23m"),c("0-5m","24-60m"),c("6-11m","12-23m"), c("6-11m","24-60m"),c("12-23m","24-60m")) #All comparisons

 
alpha_age_stratum_case<- plot_richness(PSB.sub_alpha_case_only, measures=alha_met , x="age_stratum", color="age_stratum") + 
  geom_boxplot(alpha=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = FALSE, symnum.args = symnum.args) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank()
  ) +
  geom_point(size = 2) +
  scale_color_manual(values=col_fil) +
  labs(y = "Alpha diversity") 
#ggtitle("Alpha diversity_age_stratum_cases")##
alpha_age_stratum_case



##############Beta diversity#####################################

##############Case_control_ betadiversity#####################

GP.ord <- ordinate(PSBJ.CSS.filtered, "PCoA", "bray")

beta.case_control.bray = plot_ordination(PSBJ.CSS.filtered, GP.ord, color="control_fct") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) + 
  geom_point(aes(fill = control_fct), color= "black" ,pch=21, size = 3.5) +
  scale_y_reverse() #+# Revert x-axis to match plots
  #ggtitle("case_control bray curtis") 
#ggplot2::annotate("text", x = -0.3, y = 0.3, label = "padj= 0.0001") #Add p-value from adonis

beta.case_control.bray

##Permanova _Case_control########

case_control_Permanova <- pairwise.perm.manova(phyloseq::distance(PSBJ.CSS.filtered, "bray"),PSBJ.CSS.filtered@sam_data$control_fct,p.method="fdr",nperm=9999,R2 = TRUE)

case_control_Permanova$R2.value ### R2 value
case_control_Permanova$p.value #### p-value


#### Diarrhea duration _Betadiversity ######

PSBJ.CSS.filtered_sub_bary_dia <- subset_samples(PSBJ.CSS.filtered, diarrhoea_category!="NA") ### removing  NA  #######

GP.ord <- ordinate(PSBJ.CSS.filtered_sub_bary_dia, "PCoA", "bray")

beta.diarrhoea_bary= plot_ordination(PSBJ.CSS.filtered_sub_bary_dia, GP.ord, color="diarrhoea_category", axes = 1:2) + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) + 
  geom_point(aes(fill = diarrhoea_category), color= "black" ,pch=21, size = 3.5) +
  scale_y_reverse()#+# Revert x-axis to match plots
 #ggtitle("dairroea duration bray curtis") ####
#ggplot2::annotate("text", x = -0.3, y = 0.3, label = "P < 0.001") #Add p-value from adonis

beta.diarrhoea_bary


##Permanova diarrhoea duration ################

duration_perm.manova <- pairwise.perm.manova(phyloseq::distance(PSBJ.CSS.filtered_sub_bary_dia, "bray"),PSBJ.CSS.filtered_sub_bary_dia@sam_data$diarrhoea_category,p.method="fdr",nperm=9999, R2 = TRUE)


duration_perm.manova$R2.value ## R2 value
duration_perm.manova$p.value  ### p-value

  


## Partial Constrained for diarrhea after we control age month, site, season,gender)

PSBJ.CSS.filtered_sub_rda<- subset_samples(PSBJ.CSS.filtered, diarrhoea_category!="NA")


GP.ord_dbRDA <- ordinate(PSBJ.CSS.filtered_sub_rda, "CAP","bray", formula= ~ diarrhoea_category+Condition(agemonth+site+season+gender))


p = plot_ordination(PSBJ.CSS.filtered_sub_rda, GP.ord_dbRDA, color="diarrhoea_category") +
  stat_ellipse(aes(color = diarrhoea_category),geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.line = element_line(size = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_point(aes(fill = diarrhoea_category), color= "black" ,pch=21, size = 3.5) +
  ggtitle("Constrained RDA - diarrhoea_category")

graph2ppt(p, file = "/Users/cfr965/Desktop/Getyan_script_analysis/Final_Rscriptand_figures_Deseq2_sPLS/DESeq2_cc_substudy_corrected_06_03_2022/barplot_alpha_beta_diversity",paper= "A4",width = 8, height = 6,vector.graphic = TRUE,append=TRUE)



p



###   Capscale  and ANOVA for diarrhea categories Compared  conditioned  conditioned with agemonth, site and season. 

##subset diarrhea category with no NAS##

PSBJ.CSS.filtered_sub_diarrheoa_cat<- subset_samples(PSBJ.CSS.filtered, diarrhoea_category!="NA")

##subset AD  and ProPD 

PSBJ.CSS.filtered_sub_AD_ProPD<- subset_samples(PSBJ.CSS.filtered_sub_diarrheoa_cat,diarrhoea_category!="Control")



expl_var<-c("control_fct","site","agemonth","gender","season","diarrhoea_category")



resp_var <- as.data.frame(otu_table(PSBJ.CSS.filtered_sub_AD_ProPD)) %>%
  #dplyr::select(-taxonomy) %>%
  mutate_all(as.numeric)    #### response variable the OTU Table
resp_var <- t(resp_var) ### transpose the otu table 

resp_var <- as.data.frame(resp_var) ## change to data frame.


meta_expl = as.data.frame(as.matrix(sample_data(PSBJ.CSS.filtered_sub_AD_ProPD)))

rownames <- rownames(meta_expl)  ## explanatory variable  metadata 

meta_expl <- meta_expl[rownames(resp_var),]

metadata_expl <- dplyr::select(meta_expl,expl_var)

metadata_expl


### Change the  metadata_expl to numeric and factor#


metadata_expl <- metadata_expl %>% mutate (
  control_fct = factor(control_fct), site = factor(site), agemonth = as.numeric(agemonth), gender = factor(gender),season=factor(season),diarrhoea_category=factor(diarrhoea_category))




###remove NA FROM METADATA##

metadata_expl <- metadata_expl %>% drop_na()

resp_var <- resp_var[rownames(metadata_expl),]


#Capscale AD-ProPD
capscale_AD_ProD <- capscale(resp_var ~ diarrhoea_category+ Condition(agemonth+season+gender+site) ,data =metadata_expl, dist="bray")

capscale_AD_ProD$
  
  
  Anova_capscale_AD_ProD <- anova(capscale_AD_ProD)
Anova_capscale_AD_ProD


#R2_value =SumOfSqs/ total sum of square 
#The total sum of square is model sum of square plus the residual sum of square 

R2=0.642 /(0.642 +175.877)
R2



##subset AD  and Control

PSBJ.CSS.filtered_sub_AD_Control<- subset_samples(PSBJ.CSS.filtered_sub_diarrheoa_cat,diarrhoea_category!="ProPD")



expl_var<-c("control_fct","site","agemonth","gender","season","diarrhoea_category")



resp_var <- as.data.frame(otu_table(PSBJ.CSS.filtered_sub_AD_Control)) %>%
  #dplyr::select(-taxonomy) %>%
  mutate_all(as.numeric)    #### response variable the OTU Table
resp_var <- t(resp_var) ### transpose the otu table 

resp_var <- as.data.frame(resp_var) ## change to data frame.


meta_expl = as.data.frame(as.matrix(sample_data(PSBJ.CSS.filtered_sub_AD_Control)))

rownames <- rownames(meta_expl)  ## explanatory variable  metadata 

meta_expl <- meta_expl[rownames(resp_var),]

metadata_expl <- dplyr::select(meta_expl,expl_var)
metadata_expl$diarrhoea_category


### Change the  metadata_expl to numeric and factor#


metadata_expl <- metadata_expl %>% mutate (
  control_fct = factor(control_fct), site = factor(site), agemonth = as.numeric(agemonth), gender = factor(gender),season=factor(season),diarrhoea_category=factor(diarrhoea_category))





###remove NA FROM METADATA##

metadata_expl <- metadata_expl %>% drop_na()

resp_var <- resp_var[rownames(metadata_expl),]


#Capscal AD-Control
capscale_AD_Control <- capscale(resp_var ~ diarrhoea_category+ Condition(agemonth+season+gender+site) ,data =metadata_expl, dist="bray")

capscale_AD_Control

Anova_capscale_AD_Control<- anova(capscale_AD_Control)

Anova_capscale_AD_Control

#R2_value =SumOfSqs/ total sum of square
#The total sum of square is model sum of square plus the residual sum of square 
R2=5.53/(5.43+331.69)
R2


##subset ProPD  control###



PSBJ.CSS.filtered_sub_ProPD_Control<- subset_samples(PSBJ.CSS.filtered_sub_diarrheoa_cat,diarrhoea_category!="AD")



expl_var<-c("control_fct","site","agemonth","gender","season","diarrhoea_category")



resp_var <- as.data.frame(otu_table(PSBJ.CSS.filtered_sub_ProPD_Control)) %>%
  #dplyr::select(-taxonomy) %>%
  mutate_all(as.numeric)    #### response variable the OTU Table
resp_var <- t(resp_var) ### transpose the otu table 

resp_var <- as.data.frame(resp_var) ## change to data frame.


meta_expl = as.data.frame(as.matrix(sample_data(PSBJ.CSS.filtered_sub_ProPD_Control)))

rownames <- rownames(meta_expl)  ## explanatory variable  metadata 

meta_expl <- meta_expl[rownames(resp_var),]

metadata_expl <- dplyr::select(meta_expl,expl_var)
metadata_expl$diarrhoea_category


### Change the  metadata_expl to numeric and factor#


metadata_expl <- metadata_expl %>% mutate (
  control_fct = factor(control_fct), site = factor(site), agemonth = as.numeric(agemonth), gender = factor(gender),season=factor(season),diarrhoea_category=factor(diarrhoea_category))





###remove NA FROM METADATA##

metadata_expl <- metadata_expl %>% drop_na()

resp_var <- resp_var[rownames(metadata_expl),]


#Capscal ProPD-Control
capscale_ProPD_Control <- capscale(resp_var ~ diarrhoea_category+ Condition(agemonth+season+gender+site) ,data =metadata_expl, dist="bray")

capscale_ProPD_Control

Anova_capscale_ProPD_Control<- anova(capscale_ProPD_Control)
Anova_capscale_ProPD_Control


#R2_value =SumOfSqs/ total sum of square 
#The total sum of square is model sum of square plus the residual sum of square 

R2=2.637/(2.637+194.528)
R2

####### beta diversity age_stratum #######

PSBJ.CSS.filtered_sub_bary_age_stratum <- subset_samples(PSBJ.CSS.filtered, age_stratum!="NA") ### subsetting NA#####
PSBJ.CSS.filtered_sub_bary_age_stratum@sam_data$age_stratum<- factor(PSBJ.CSS.filtered_sub_bary_age_stratum@sam_data$age_stratum,levels=c("0-5m", "6-11m","12-23m", "24-60m")) 

GP.ord <- ordinate(PSBJ.CSS.filtered_sub_bary_age_stratum, "PCoA", "bray")

beta.age_stratum_bary= plot_ordination(PSBJ.CSS.filtered_sub_bary_age_stratum, GP.ord, color="age_stratum") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) + 
  geom_point(aes(fill = age_stratum), color= "black" ,pch=21, size = 3.5) +
  scale_y_reverse()#+# Revert x-axis to match plots
  #ggtitle(" age_stratum bray curtis") ####
#ggplot2::annotate("text", x = -0.3, y = 0.3, label = "P < 0.001") #Add p-value from adonis

beta.age_stratum_bary



age_stratum_pairwise.perm.manova<- pairwise.perm.manova(phyloseq::distance(PSBJ.CSS.filtered_sub_bary_age_stratum, "bray"),PSBJ.CSS.filtered_sub_bary_age_stratum@sam_data$age_stratum,p.method="fdr",nperm=9999,R2 = TRUE)

 
age_stratum_pairwise.perm.manova$R2.value ### R2 value 
age_stratum_pairwise.perm.manova$p.value  ### p-value


#### beta diversity in cases and controls  strartfied  by age stratum####  
             
## subsetting Controls ########

PSB.sub_beta_controls_only<- subset_samples(PSBJ.filtered,control_fct!= "Case") ### controls 


PSB.sub_beta_controls_only@sam_data$age_stratum<- factor(PSB.sub_beta_controls_only@sam_data$age_stratum,levels=c("0-5m", "6-11m","12-23m", "24-60m")) 

GP.ord <- ordinate(PSB.sub_beta_controls_only, "PCoA", "bray")

beta.age_stratum_bary_controls= plot_ordination(PSB.sub_beta_controls_only, GP.ord, color="age_stratum") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) + 
  geom_point(aes(fill = age_stratum), color= "black" ,pch=21, size = 3.5) +
  scale_x_reverse()#+# Revert x-axis to match plots
#ggtitle(" age_stratum bray curtis") ####
#ggplot2::annotate("text", x = -0.3, y = 0.3, label = "P < 0.001") #Add p-value from adonis

beta.age_stratum_bary_controls




##Permanova age stratum controls ########

control_age_stratum_pairwise.perm.manova<- pairwise.perm.manova(phyloseq::distance(PSB.sub_beta_controls_only, "bray"),PSB.sub_beta_controls_only@sam_data$age_stratum,p.method="fdr",nperm=9999,R2 = TRUE)


control_age_stratum_pairwise.perm.manova$R2.value ### R2 value
control_age_stratum_pairwise.perm.manova$p.value ### p-value
  
## subset cases  beta diversity stratified by age strartum  ##

PSB.sub_beta_cases_only<- subset_samples(PSBJ.filtered,control_fct!= "Control") ### cases 


PSB.sub_beta_cases_only@sam_data$age_stratum<- factor(PSB.sub_beta_cases_only@sam_data$age_stratum,levels=c("0-5m", "6-11m","12-23m", "24-60m")) 

GP.ord <- ordinate(PSB.sub_beta_cases_only, "PCoA", "bray")

beta.age_stratum_bary_cases= plot_ordination(PSB.sub_beta_cases_only, GP.ord, color="age_stratum") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, size = 1) +
  scale_fill_manual(values= col_fil) +
  scale_color_manual(values = col_fil) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) + 
  geom_point(aes(fill = age_stratum), color= "black" ,pch=21, size = 3.5) +
  scale_x_reverse()#+# Revert x-axis to match plots
#ggtitle(" age_stratum bray curtis") ####
#ggplot2::annotate("text", x = -0.3, y = 0.3, label = "P < 0.001") #Add p-value from adonis

beta.age_stratum_bary_cases



##Permanova age stratum cases only  ####

case_age_stratum_pairwise.perm.manova<- pairwise.perm.manova(phyloseq::distance(PSB.sub_beta_cases_only, "bray"),PSB.sub_beta_cases_only@sam_data$age_stratum,p.method="fdr",nperm=9999,R2 = TRUE)

case_age_stratum_pairwise.perm.manova$p.value
case_age_stratum_pairwise.perm.manova$R2.value




#####  Effect size of demographic an clinical variables on the Gut microbiome  ###

###Distance based redundant analysis(db-RDA) for the entire cohort ( cases and controls ######


sample_data(PSBJ.CSS.filtered)$histad<-as.factor(sample_data(PSBJ.CSS.filtered)$histad)## to change the numeric variable to factor

sample_data(PSBJ.CSS.filtered)$hospitaladmission<-revalue(sample_data(PSBJ.CSS.filtered)$histad,c("0"="no","1"="yes","3"="yes")) ## reclassifying the histad
sample_data(PSBJ.CSS.filtered)$malnutrition<-revalue(sample_data(PSBJ.CSS.filtered)$am,c("SAM"="Yes","MAM"="Yes","Normal"="No")) ## reclassifying the am 


expl_var <-  c("control_fct","site","malnutrition","agemonth","gender","muac_fct","muac","stunting","biloed" ,"caretaker_fct","crowding","crowding_dich","matedu_16","matedu_trip","asset_index","animal_any","histab","histab_type","histdiar","sanitation","improvedsan_fct","water_trip","watertre_fct","histdiar","histvisi_dich", "weaned_early","prematur_fct","modedel_fct","histad","hospitaladmission","breastmilk_now","histmalr","season")  ### create a vector of clinical variables 



 
resp_var <- as.data.frame(otu_table(PSBJ.CSS.filtered)) %>%
  #dplyr::select(-taxonomy) %>%
  mutate_all(as.numeric)    #### response variable the OTU Table
resp_var <- t(resp_var) ### transpose the otu table 

resp_var <- as.data.frame(resp_var) ## change to data frame.


meta_expl = as.data.frame(as.matrix(sample_data(PSBJ.CSS.filtered)))

rownames <- rownames(meta_expl)  ## explanatory variable  metadata 

meta_expl <- meta_expl[rownames(resp_var),]

metadata_expl <- dplyr::select(meta_expl,expl_var) ## selecting variables 

metadata_expl   ### this is the clinical variables used in the db-RDA


### Change the  metadata_expl to numeric and factor accordingly #######


metadata_expl <- metadata_expl %>% mutate (
  control_fct = factor(control_fct), site = factor(site), agemonth = as.numeric(agemonth), gender = factor(gender),muac= as.numeric(muac), muac_fct = factor(muac_fct),  caretaker_fct = factor(caretaker_fct), malnutrition=factor(malnutrition), crowding = factor(crowding), crowding_dich=factor(crowding_dich), matedu_16 = as.numeric(matedu_16), matedu_trip=factor(matedu_trip), asset_index = as.numeric(asset_index), animal_any=factor(animal_any),histab = factor(histab), histab_type = factor(histab_type), sanitation=factor(sanitation), improvedsan_fct=factor(improvedsan_fct), water_trip=factor(water_trip), watertre_fct=factor(watertre_fct), histdiar = factor(histdiar), histvisi_dich = factor(histvisi_dich), weaned_early = factor(weaned_early), prematur_fct = factor(prematur_fct), modedel_fct = factor(modedel_fct) ,histad =factor(histad), hospitaladmission=factor(hospitaladmission), breastmilk_now = factor(breastmilk_now),histmalr = factor(histmalr), season=factor(season))


### create new  other variables from the metadata ###


#### Create WAM-index variable from (water/sanitation,Asset,Maternal education)
metadata_expl <- metadata_expl%>%
  mutate(sanitation_numeric = case_when (improvedsan_fct =="improved" ~ 4,improvedsan_fct == "unimproved" ~ 0)) ## assign numeric value for for sanitation 


metadata_expl<- metadata_expl%>%
  mutate(water_numeric = case_when (water_trip =="public tap" ~ 4, water_trip == "private tap" ~ 4,water_trip == "Unimproved/borehole/p.spring" ~ 0)) ## assign numeric value for for sanitation 


metadata_expl$water_sanitation <- metadata_expl$sanitation_numeric+ metadata_expl$water_numeric ###adding the values and create new variable called water_sanitation

metadata_expl$maternal_education <- (metadata_expl$matedu_16)/2


metadata_expl$WAM <- (metadata_expl$asset_index+metadata_expl$water_sanitation+metadata_expl$maternal_education)/32##  create WAM variable 

metadata_expl <- metadata_expl %>%
  mutate(WAM_index = case_when (WAM < 0.5  ~ "under_0.5", WAM >= 0.5  ~ "overeq_0.5")) #### WAM_index category 



##Drop variables that are not used db-RDA ####

metadata_expl <- select (metadata_expl,-c(muac_fct,muac,crowding,matedu_16,asset_index,improvedsan_fct,watertre_fct,water_numeric, maternal_education,matedu_trip,histab_type,sanitation, water_trip,sanitation_numeric,water_sanitation,histad,WAM))


##### Run db-RDA  model for thh whole cohort #### 

metadata_expl <- metadata_expl %>% drop_na() ## drop NA from metadta

resp_var <- resp_var[rownames(metadata_expl),] ## response variabler


db_RDA <- dbrda(resp_var ~ control_fct+agemonth+malnutrition+gender+stunting+biloed+caretaker_fct+ crowding_dich+animal_any+histab+histdiar+histvisi_dich+weaned_early+prematur_fct+modedel_fct+hospitaladmission +breastmilk_now+ histmalr+WAM_index+ Condition(site+season), data = metadata_expl,dist= "bray", na.action= na.omit) ### run the model

db_RDA

anova_db_RDA <- anova(db_RDA)  ### anova of the model

anova_db_RDA

anova_db_RDA_terms <-anova(db_RDA,by="term") ## run the anova of  individual variables in the model
anova_db_RDA_terms


## ########prepare for plotting #####

anova_db_RDA_terms_data <- as.data.frame(anova_db_RDA_terms)## change to dataframe



#Add adjusted p-values
anova_db_RDA_terms_data$`Pr(>F).adj` <- p.adjust(anova_db_RDA_terms_data$`Pr(>F)`, method = 'BH', n = ncol (metadata_expl))


#######Calculate R2  value form the  SumOfSqs  and add to the  the dataframe.

anova_db_RDA_terms_data$SumOfSqs_plus_residual <- anova_db_RDA_terms_data$SumOfSqs+195.213  

anova_db_RDA_terms_data$R2 <- anova_db_RDA_terms_data$SumOfSqs/anova_db_RDA_terms_data$SumOfSqs_plus_residual ## Calculate R2



#Plotting ordered R2
anova_db_RDA_terms_data <- as.data.frame(anova_db_RDA_terms_data) %>%
  tibble::rownames_to_column("Variable")  %>%
  filter(!Variable %in% c("Residual")) %>%
  arrange(R2)

anova_db_RDA_terms_data.melt <- tidyr::pivot_longer(anova_db_RDA_terms_data,cols=c("R2"), names_to='measure', values_to="value") %>%
  add_significance(p.col = "Pr(>F).adj", output.col = "PrF.sym",
                   cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                   symbols = c("****", "***", "**", "*", "")
  )



fancy.dbRDA.bar.all_cohort <- ggbarplot(anova_db_RDA_terms_data.melt,
                                        x="Variable",
                                        y="value",
                                        #color = "Variable",
                                        orientation = c("horizontal"),
                                        color = NA,
                                        fill = "measure",
                                        position=position_dodge(0.7),
                                        label = anova_db_RDA_terms_data.melt$PrF.sym,
                                        lab.size = 9,
                                        lab.vjust = 0.7,
                                        lab.hjust = -0.3
                                        
) +
  
  ggtitle("drRDA dissimilarity effect size_all_cohort") +
  ylab(expression(R^2~"value")) +
  scale_fill_jco()

fancy.dbRDA.bar.all_cohort



write.csv(anova_db_RDA_terms_data.melt,"drDRDA_all_cohort.csv") ## exported as csv



###db-RDA Diarrheoa cases ony  ###


expl_var <-  c("site","malnutrition","agemonth","gender","muac_fct","diarrhoea_category","stool", "curblo","histors","muac","stunting","biloed" ,"caretaker_fct","crowding","crowding_dich","matedu_16","matedu_trip","asset_index","animal_any","histab","histab_type","histdiar","sanitation","improvedsan_fct","water_trip","watertre_fct","histdiar","histvisi_dich", "weaned_early","prematur_fct","modedel_fct","histad","hospitaladmission","breastmilk_now","histmalr","season") ## create vector for variable diarrhea cases only 


#### dr-RDA  diarrhea cases only ######

PSBJ.CSS.filtered <- subset_samples(PSBJ.CSS.filtered,diarrhoea_category!="Control")### remove all controls and keep only diarrhea cases. 

resp_var <- as.data.frame(otu_table(PSBJ.CSS.filtered)) %>%
  #dplyr::select(-taxonomy) %>%
  mutate_all(as.numeric)    #### response variable the OTU Table
resp_var <- t(resp_var) ### transpose the otu table 

resp_var <- as.data.frame(resp_var) ## change to data frame.


meta_expl = as.data.frame(as.matrix(sample_data(PSBJ.CSS.filtered)))

rownames <- rownames(meta_expl)  ## explanatory variable  metadata 

meta_expl <- meta_expl[rownames(resp_var),]

metadata_expl <- dplyr::select(meta_expl,expl_var)

metadata_expl


### Change the  metadata_expl to numeric and factor#

metadata_expl <- metadata_expl %>% mutate (
   site = factor(site), agemonth = as.numeric(agemonth), gender = factor(gender),muac= as.numeric(muac), muac_fct = factor(muac_fct),  caretaker_fct = factor(caretaker_fct), malnutrition=factor(malnutrition), crowding = factor(crowding), crowding_dich=factor(crowding_dich), matedu_16 = as.numeric(matedu_16), matedu_trip=factor(matedu_trip), asset_index = as.numeric(asset_index), animal_any=factor(animal_any),histab = factor(histab), histab_type = factor(histab_type), sanitation=factor(sanitation), improvedsan_fct=factor(improvedsan_fct), water_trip=factor(water_trip), watertre_fct=factor(watertre_fct), histdiar = factor(histdiar), histvisi_dich = factor(histvisi_dich), weaned_early = factor(weaned_early), prematur_fct = factor(prematur_fct), modedel_fct = factor(modedel_fct) ,histad =factor(histad), hospitaladmission=factor(hospitaladmission), breastmilk_now = factor(breastmilk_now),histmalr = factor(histmalr), season=factor(season), diarrhoea_category= factor(diarrhoea_category),stool= factor(stool), curblo=factor(curblo),histors= factor(histors))



## create new variable MUAC_undereq125_and_age_above6months

metadata_expl<- metadata_expl%>%
  mutate(muac_undereq125 = case_when (muac <= 125 & agemonth > 6 ~ "undereq_125", muac > 125 & agemonth > 6 ~ "over_125"))


##Measuring socioeconomic status in multidisciplinary: results from the eight-country MAL-ED study us these paper to produce WAM (water/sanitation,Asset,Maternal education) variable blow. 


metadata_expl <- metadata_expl%>%
  mutate(sanitation_numeric = case_when (improvedsan_fct =="improved" ~ 4,improvedsan_fct == "unimproved" ~ 0)) ## assign numeric value for for sanitation 


metadata_expl<- metadata_expl%>%
  mutate(water_numeric = case_when (water_trip =="public tap" ~ 4, water_trip == "private tap" ~ 4,water_trip == "Unimproved/borehole/p.spring" ~ 0)) ## assign numeric value for for sanitation 


metadata_expl$water_sanitation <- metadata_expl$sanitation_numeric+ metadata_expl$water_numeric ###adding the values and create new variable called water_sanitation

metadata_expl$maternal_education <- (metadata_expl$matedu_16)/2


metadata_expl$WAM <- (metadata_expl$asset_index+metadata_expl$water_sanitation+metadata_expl$maternal_education)/32##  creat WAM variable 

metadata_expl <- metadata_expl %>%
  mutate(WAM_index = case_when (WAM < 0.5  ~ "under_0.5", WAM >= 0.5  ~ "overeq_0.5"))




##drop variables that are not important####

metadata_expl <- select (metadata_expl,-c(muac_fct,muac,crowding,matedu_16,asset_index,improvedsan_fct,watertre_fct,water_numeric, maternal_education,matedu_trip,histab_type,sanitation, water_trip,sanitation_numeric,water_sanitation,histad,WAM))


#####db-DRA diarrhea cases only #########


metadata_expl <- metadata_expl %>% drop_na()

resp_var <- resp_var[rownames(metadata_expl),]


db_RDA_case <- dbrda(resp_var ~agemonth+malnutrition+gender+stunting+biloed+caretaker_fct+ crowding_dich+animal_any+histab+histdiar+histvisi_dich+weaned_early+prematur_fct+modedel_fct+hospitaladmission +breastmilk_now+ histmalr+WAM_index+diarrhoea_category+stool+curblo+histors + Condition(site+season), data = metadata_expl,dist= "bray", na.action= na.omit)    ### the full model
db_RDA_case



anova_db_RDA_case<- anova(db_RDA_case) ### 

anova_db_RDA_case



anova_db_RDA_case_terms<- anova(db_RDA_case,by="term") ## anova for individual all variables
anova_db_RDA_case_terms


##prepare for ploting 

anova_db_RDA_case_terms_data <- as.data.frame(anova_db_RDA_case_terms)## change it  to data frame


#Add adjusted p-values
anova_db_RDA_case_terms_data$`Pr(>F).adj` <- p.adjust(anova_db_RDA_case_terms_data$`Pr(>F)`, method = 'BH', n = ncol (metadata_expl))


### calculate R2  value form the  SumOfSqs  and add to  the data frame.

anova_db_RDA_case_terms_data$SumOfSqs_plus_residual <- anova_db_RDA_case_terms_data$SumOfSqs+115.761 

anova_db_RDA_case_terms_data$R2 <- anova_db_RDA_case_terms_data$SumOfSqs/anova_db_RDA_case_terms_data$SumOfSqs_plus_residual ## Calculate R2

#Show table
#knitr::kable(anova.dat)

#Plotting ordered R2
anova_db_RDA_case_terms_data <- as.data.frame(anova_db_RDA_case_terms_data) %>%
  tibble::rownames_to_column("Variable")  %>%
  filter(!Variable %in% c("Residual")) %>%
  arrange(R2)

anova_db_RDA_case_terms_data.melt <- tidyr::pivot_longer(anova_db_RDA_case_terms_data,cols=c("R2"), names_to='measure', values_to="value") %>%
  add_significance(p.col = "Pr(>F).adj", output.col = "PrF.sym",
                   cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                   symbols = c("****", "***", "**", "*", "")
  )



fancy.dbRDA.bar_diarrhea_only <- ggbarplot(anova_db_RDA_case_terms_data.melt,
                                           x="Variable",
                                           y="value",
                                           #color = "Variable",
                                           orientation = c("horizontal"),
                                           color = NA,
                                           fill = "measure",
                                           position=position_dodge(0.7),
                                           label = anova_db_RDA_case_terms_data.melt$PrF.sym,
                                           lab.size = 9,
                                           lab.vjust = 0.7,
                                           lab.hjust = -0.3
                                           
) +
  #stat_pvalue_manual(
  #  adonis.melt,  label = "PrF.sym", tip.length = 0.01
  #) +
  #geom_hline(yintercept = 0.05) +
  ggtitle("drRDA decomposed dissimilarity effect size: non-collinear variables diarrhoea cases only") +
  ylab(expression(R^2~"value")) +
  scale_fill_jco()

fancy.dbRDA.bar_diarrhea_only



write.csv(anova_db_RDA_case_terms_data.melt,"drDRDA_Diarrhoeacaseonly.csv") ## exported as csv



### Logistic regression analysis of diarrhea ( Yes or No) for demographic and clinical variables 


sample_data(PSBJ.filtered)$case_control_cat<- revalue(sample_data(PSBJ.filtered)$control_fct,c("Case"= "1", "Control"="0"))  ##### case and control


## Run the logistic regeneration model  model case_control


table(datalog$case_control_cat) ##check the data 

datalog### metadata

final_model <- glm(case_control_cat ~ caretaker_fct+ modedel_fct+prematur_fct+breastmilk_now+crowding_dich+muac_undereq125+WAM_index+weaned_early+histab+histdiar+histmalr+histvisi_dich+hospitaladmission+animal_any+agemonth+site+season+gender, family = binomial(link = "logit"), data = datalog) ### include clinical important variables 

summary(final_model)
mod <- exp(coef(final_model))  ## to the adds ratio
mod
exp(confint(final_model))  ## confidence interval 



### Logistic regression prolonged or persistent diarrhea (ProPD)  and Acute diarrhea (AD)


PS_CASE_only<- subset_samples(PSBJ.filtered,control_fct!="Control") ## subset diarrhea cases only 
PS_CASE_only

sample_data(PS_CASE_only)$diarrhoea_AD_ProPD<- revalue(sample_data(PS_CASE_only)$diarrhoea_category,c( "ProPD"="1", "AD"= "0")) ##  


## Final Model logist regression ProPD_AD ###

diarrhoea_AD_ProPD ### metadata

Final_AD_ProPD <- glm(diarrhoea_AD_ProPD ~stunting+histors+histab+histdiar+muac_undereq125+WAM_index+weaned_early +prematur_fct+hospitaladmission+histmalr+histvisi_dich+caretaker_fct+breastmilk_now+animal_any+stool_category+crowding_dich+curblo+modedel_fct+agemonth+site+season+gender, family = binomial(link = "logit"), data = dataProPD)  ## include clinical important variables  

summary(Final_AD_ProPD) 

exp(coef(Final_AD_ProPD))  ## to the adds ratio
exp(confint(Final_AD_ProPD)) ## confidence interval 

