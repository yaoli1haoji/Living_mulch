if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") 
#install.packages ("rlang")
#remove.packages("data.table")
#remove.packages("ggplot2")
#install.packages ("data.table")
#install.packages ("ggplot2")
install_github("vqv/ggbiplot")
library(qiime2R)
library(ggplot2)
library(tidyverse)
library(devtools)
library(gridExtra)
library(ggpubr)
library(rstatix)

#set the working directory 
setwd("/home/hl46161/publish_living_mulch/")

#read in metadata and shannon index 
metadata <- read.table("living_mulch_data_metafile_massalin.tsv", sep="\t",header = TRUE)
shannon<-read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") 
shannon

##read in evenness artifact
eveness <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/evenness_vector.qza")
eveness <-eveness$data %>% rownames_to_column("SampleID") 

##read in observed otu artifact
observed_otu <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/observed_features_vector.qza")
observed_otu <- observed_otu$data %>% rownames_to_column("SampleID")
#merge the data 
metadata <- merge(shannon,metadata, by = "SampleID")
metadata <- merge(eveness,metadata, by = "SampleID")
metadata <- merge(observed_otu,metadata, by = "SampleID")
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
shannon$data 

#################

# install.packages("rJava")
#install.packages("UpSetR")
#install.packages("tidyverse")
#install.packages("venneuler")
#install.packages("grid")

library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)


############plot the pielou evenness of each treatment seperate by sampling date 
gg <-ggplot(metadata,aes(x=Treatment, y=pielou_evenness)) + 
  facet_grid(~`Date_taken`) +
  xlab("Treatment") + 
  facet_grid(~`Date_taken`) +
  ylab("Pielou_Evenness") + 
  geom_point(aes(color=Treatment),size=4) +
  theme(
    legend.text = element_text(color = "black", size = 20),   #set the text size to be 20 and bold the text
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20, face="bold")
  ) + 
  theme(
    panel.grid.major = element_blank(),   ##remove grid line and background 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 20)) +
  geom_line() 

gg  

#save the plot in /home/hl46161/publish_living_mulch/ directroy 
ggsave("Pielou_Evenness.pdf", height=8, width=8, device="pdf")
ggsave("Pielou_Evenness.png", height=8, width=8, device="png")

############plot the shannon diversity of each treatment seperate by sampling date 
gg <-ggplot(metadata,aes(x=Treatment, y=shannon_entropy)) + 
  facet_grid(~`Date_taken`) +
  xlab("Treatment") + 
  facet_grid(~`Date_taken`) +
  ylab("shannon_entropy") + 
  geom_point(aes(color=Treatment),size=4) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20, face="bold")
  ) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 20)) +
  geom_line() 

gg  

#save the plot in /home/hl46161/publish_living_mulch/ directroy 
ggsave("shannon_entropy.pdf", height=8, width=8, device="pdf")
ggsave("shannon_entropy.png", height=8, width=8, device="png")

################################ Try to use linear regression to determine the effect of the treatment and sampling date on alpha diversity 

treatment_order <- metadata$Treatment
#change the reference level of treatment
treatment_order <- relevel(metadata$Treatment, ref="NoCover")

#run linear regression of shannon index on treatment and date taken and their interaction
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order + Date_taken, data = metadata)

#summary the model 
summary(lm_treatment_date_relation_shannon)
### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
plot(lm_treatment_date_relation_shannon)


#######################################################################
###assumption of model above is not satisfied, drop 11th and 12th samples according to plots of previous function
### 11th and 12th sample are 
treatment_order <- treatment_order[-c(11,12)]

#### rerun the model after filteration 
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order*Date_taken, data = metadata[-c(11,12),])

######summary the model
summary(lm_treatment_date_relation_shannon)
anova(lm_treatment_date_relation_shannon)
### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
plot(lm_treatment_date_relation_shannon)

##since the iteraction is not significant, we can test treatment and Date_taken sperately using type 2 anova from car R Package
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order  + Date_taken, data = metadata[-c(11,12),])

######summary the model 
summary(lm_treatment_date_relation_shannon)
Anova(lm_treatment_date_relation_shannon,type=2)
### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
plot(lm_treatment_date_relation_shannon)

################################### run the model of evenness on treatment and date_taken, and their interaction
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_evenness)
##########significant differences exist in iteraction of living mulch and June 28th and August 31 
### change the model to run type three Anova
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order*Date_taken, data = metadata[-c(11,12),],contrasts=list(treatment_order=contr.sum, Date_taken=contr.sum))
summary(lm_treatment_date_relation_evenness)
Anova(lm_treatment_date_relation_evenness,type=3)
gvlma::gvlma(lm_treatment_date_relation_evenness)
#########
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_evenness)
Anova(lm_treatment_date_relation_shannon,type=2)
gvlma::gvlma(lm_treatment_date_relation_evenness)
plot(lm_treatment_date_relation_evenness)


#######################################################
###no significant effect of either treatment and sampling date on # of ASV observed 
lm_treatment_date_relation_observed_otu <- lm(formula = observed_features ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_observed_otu)


################################################################33
library(emmeans)
library(multcomp)
########## run a linear regression of soil physical characteristic on treatment and date_taken 
### The goal is to see whether treatment affect the soil physical characteritics 

lm_pH <- lm(formula = pH ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_pH) # no significant iteraction, so test treatment and date seperately 

lm_pH <- lm(formula = pH ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_pH)

##### no significant differnce between treatment
emm1 = emmeans(lm_pH, specs = pairwise ~ treatment_order)
emm1$contrasts


#######################################################

lm_LBCEQ <- lm(formula = LBCEQ ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_LBCEQ) # significant iteraction exist between CrimsonClove and 2018-06-28

lm_LBCEQ <- lm(formula = LBCEQ ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_LBCEQ)

######## All three cover crop treatment significant different from nocover 
emm1 = emmeans(lm_LBCEQ, specs = pairwise ~ treatment_order)
emm1$contrasts

## generate compact letter display 
cld <- cld(glht(lm_LBCEQ, linfct=mcp(treatment_order = "Tukey")))
cld

#find out the max value of each treatment of LBCEQ
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(LBCEQ))
value_max
### add the compact letter display to max value and organize into data frame 
value_max$letters <- c("A","A","A","B")
value_max <- as.data.frame(value_max)
###### create a treatment column that match with metadata so that the value max dataframe can match and be plotted on the graph
rownames(value_max) <- value_max$Treatment
value_max

#plot the LBCEQ value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=LBCEQ)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 10+max_value,vjust=0),label=value_max$letters) # input the compact letter display on the top of max value

ggsave("LBCEQ.pdf", height=8, width=10, device="pdf")
ggsave("LBCEQ.png", height=8, width=10, device="png")

##################################################################################
### no significant difference of treatment in CEC  
lm_CEC <- lm(formula = CEC ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_CEC)
lm_CEC <- lm(formula = CEC ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_CEC)

############################ 
emm1 = emmeans(lm_CEC, specs = pairwise ~ treatment_order)
emm1$contrasts

#######interesting 

#################
####the model is not significant 
###no need to further explore  
lm_Ca <- lm(formula = Ca ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Ca)
lm_Ca <- lm(formula = Ca ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Ca)

emm1 = emmeans(lm_Ca, specs = pairwise ~ treatment_order)
emm1$contrasts

####################
lm_K <- lm(formula = K ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_K)
####### no significant iteraction between treatment and date taken
lm_K <- lm(formula = K ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_K)

######## All three cover crop treatment significant different from nocover 
emm1 = emmeans(lm_K, specs = pairwise ~ treatment_order)
emm1$contrasts


cld(glht(lm_K, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(K))
value_max
value_max$letters <- c("AB","A","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the K value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=K)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 1.05*max_value,vjust=0),label=value_max$letters)


ggsave("K.pdf", height=8, width=10, device="pdf")
ggsave("K.png", height=8, width=10, device="png")


################

lm_Mg <- lm(formula = Mg ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Mg)  ### significant iteraction exist between LivingMulch and 2018-06-28

lm_Mg <- lm(formula = Mg ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Mg) ### Living mulch and Crimsonclover significantly different from Nocover 

emm1 = emmeans(lm_Mg, specs = pairwise ~ treatment_order)
emm1$contrasts

## generate compact letter display 
cld(glht(lm_Mg, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Mg))
value_max$letters <- c("A","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Mg value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=Mg)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 5 + max_value,vjust=0),label=value_max$letters)

ggsave("Mg.pdf", height=8, width=10, device="pdf")
ggsave("Mg.png", height=8, width=10, device="png")


###################################

lm_Zn <- lm(formula = Zn ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Zn)  ### no significant iteraction exist 

lm_Zn <- lm(formula = Zn ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Zn)  ######## all three cover crop significantly different 

emm1 = emmeans(lm_Zn, specs = pairwise ~ treatment_order)
emm1$contrasts

## generate compact letter display 
cld(glht(lm_Zn, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Zn))
value_max$letters <- c("B","B","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Zn value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=Zn)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 0.2 + max_value,vjust=0),label=value_max$letters)

ggsave("Zn.pdf", height=8, width=10, device="pdf")
ggsave("Zn.png", height=8, width=10, device="png")


####################################

lm_Amon <- lm(formula = Amon ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Amon)  ### no significant iteraction exist 

lm_Amon <- lm(formula = Amon ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Amon)  #CereaRye amd Living mulch significant different 

emm1 = emmeans(lm_Amon, specs = pairwise ~ treatment_order)
emm1$contrasts
## generate compact letter display 
cld(glht(lm_Amon, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of Amon
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Amon))
value_max$letters <- c("AB","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Amon value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=Amon)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 1 + max_value,vjust=0),label=value_max$letters)

ggsave("Amon.pdf", height=8, width=10, device="pdf")
ggsave("Amon.png", height=8, width=10, device="png")


###########################
lm_Nit <- lm(formula = Nit ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Nit) ## no significant correlatiob between treatment and sampling date 
lm_Nit <- lm(formula = Nit ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Nit) ### no tratment is significantly different, but sampling date are significantly different, which is understandable  

emm1 = emmeans(lm_Nit, specs = pairwise ~ treatment_order)
emm1$contrasts

##############################

lm_N <- lm(formula = N ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_N)  #no significant iteraction exist 

lm_N <- lm(formula = N ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_N) # all three treatment are significantly different from no cover 

emm1 = emmeans(lm_N, specs = pairwise ~ treatment_order)
emm1$contrasts
## generate compact letter display 
cld(glht(lm_N, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of N
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(N))
value_max$letters <- c("B","BC","C","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the N value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=N)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y =  1.05*max_value,vjust=0),label=value_max$letters)

ggsave("N.pdf", height=8, width=10, device="pdf")
ggsave("N.png", height=8, width=10, device="png")



#################################

lm_TOC <- lm(formula = TOC ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_TOC) # no significant iteraction

lm_TOC <- lm(formula = TOC ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_TOC) #Living mulch is significant on TOC

emm1 = emmeans(lm_TOC, specs = pairwise ~ treatment_order)
emm1$contrasts

## generate compact letter display 
cld(glht(lm_TOC, linfct=mcp(treatment_order = "Tukey")))
#find out the max value of each treatment of TOC
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(TOC))
value_max$letters <- c("A","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the TOC value based on treatment and colored by date 
ggplot(metadata,aes(x=treatment,y=TOC)) +  geom_line() +  geom_point(aes(col=Date_taken),size=6) +  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y =  1.05*max_value,vjust=0),label=value_max$letters)

ggsave("TOC.pdf", height=8, width=10, device="pdf")
ggsave("TOC.png", height=8, width=10, device="png")

###############################################################

#Test the association between soil physical data and alpha diversity, first start with evenness

metadata_late = subset(metadata[-c(11,12),],Date_taken == "2018-08-31")

n <- 20
####### for loop the soil physical data column with evenness entropy 
lm_evenness_soil_august <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_late[,x], data = metadata_late))
######### summarize the model
evenness_soil_august_summaries <- lapply(lm_evenness_soil_august, summary)
evenness_soil_august_summaries
########## The 2 3 5 7 13 column show significant correlation 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
##############Therefore pH BS Ca Mg N significantly correlates with evenness 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_august, gvlma::gvlma)
assumption_check
########### 

####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_august <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_late[,x], data = metadata_late))
######### summarize the model
shannon_soil_august_summaries <- lapply(lm_shannon_soil_august, summary)
shannon_soil_august_summaries
## pH is significant 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_august, gvlma::gvlma)
########## For august no soil physical data is associated with shannon entropy 

#######################################################################

metadata_mid = subset(metadata[-c(11,12),],Date_taken == "2018-06-28")
n <- 20
####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_mid <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_mid[,x], data = metadata_mid))
######### summarize the model
shannon_soil_mid_summaries <- lapply(lm_shannon_soil_mid, summary)
shannon_soil_mid_summaries
######### 10 13 models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"
### therefore nitrogen and ammonia is signifiant 
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
############## N and AMon is strongly associated with shannon entropy in June 
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_mid, gvlma::gvlma)
assumption_check
########## 

####### for loop the soil physical data column with pielou evenness 
lm_evenness_soil_mid <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_mid[,x], data = metadata_mid))
######### summarize the model 
evenness_soil_mid_summaries <- lapply(lm_evenness_soil_mid, summary)
evenness_soil_mid_summaries
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
#"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
########## N and AMon is strongly associated with pielou evenness in June 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_mid, gvlma::gvlma)
assumption_check

###############################################

metadata_early = subset(metadata[-c(11,12),],Date_taken == "2018-05-21")
n <- 20
####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_early <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_early[,x], data = metadata_early))
lm_shannon_soil_early <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_early[,x], data = metadata_early))
######### summarize the model
shannon_soil_early_summaries <- lapply(lm_shannon_soil_early, summary)
shannon_soil_early_summaries
######### 9th models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N" 
#### therefore Zn is significant in shannon
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_early, gvlma::gvlma)
assumption_check
########## 





###############################################################
##plot mantel test matrix 

calc_distance_vector <- function(DM){
  
  rownames(DM) <- DM$X
  ##delete previous sample names column
  DM$X <-NULL
  ## reorder the dataframe by rownames and column names 
  DM <-DM[order(row.names(DM)),  order(colnames(DM))]
  
  #set up loop numbers 
  row_loop=nrow(DM)
  print(row_loop)
  column_loop=ncol(DM)-1
  print(column_loop)
  ##set up an emptry vector
  unifrac_distance_vector <- c(0)
  ###

  
  ###### loop trhough the lower triangle of the matrix and add elements to the vector 
  while(row_loop <=nrow(DM) & row_loop >1){
    for(i in column_loop:1){
      #print(i)
      #print(DM[row_loop,i])
      unifrac_distance_vector <- append(unifrac_distance_vector,DM[row_loop,i])
    }
    column_loop=column_loop-1
    row_loop = row_loop-1
  }
  ### remove first zero 
  unifrac_distance_vector <- unifrac_distance_vector[2:length(unifrac_distance_vector)]
  unifrac_distance_vector
  return(unifrac_distance_vector)
}
  
######################################################################3
## read in august weighted unifrac distance 
August_weighted_unifrac_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-weighted_unifrac_distance-matrix.tsv",sep="\t",header=TRUE)
August_weighted_unifrac_distance_vector <- calc_distance_vector(August_weighted_unifrac_distance)
August_weighted_unifrac_distance_vector

August_unweighted_unifrac_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-unweighted_unifrac_distance-matrix.tsv",sep="\t",header=TRUE)
August_unweighted_unifrac_distance_vector <- calc_distance_vector(August_unweighted_unifrac_distance)
August_unweighted_unifrac_distance_vector

August_bray_curtis_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-bray-curtis_distance.tsv",sep="\t",header=TRUE)
August_bray_curtis_distance_vector <- calc_distance_vector(August_weighted_unifrac_distance)
August_bray_curtis_distance_vector

#####
## read in august N unifrac distance 
August_N_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-N-distance-matrix.tsv",sep="\t",header=TRUE)
August_N_distance_vector <- calc_distance_vector(August_N_distance)
August_N_distance_vector

##############################
## read in august N unifrac distance 
August_LBCEQ_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-LBCEQ-distance-matrix.tsv",sep="\t",header=TRUE)
August_LBCEQ_distance_vector <- calc_distance_vector(August_LBCEQ_distance)
August_LBCEQ_distance_vector

########################################

August_data <- data.frame(August_N_distance_vector,August_LBCEQ_distance_vector,August_weighted_unifrac_distance_vector,August_unweighted_unifrac_distance_vector,August_bray_curtis_distance_vector)

summary(lm(August_weighted_unifrac_distance_vector ~ August_N_distance_vector,data=August_data))
summary(lm(August_unweighted_unifrac_distance_vector ~ August_N_distance_vector,data=August_data))
summary(lm(August_bray_curtis_distance_vector ~ August_N_distance_vector,August_data))


summary(lm(August_weighted_unifrac_distance_vector ~ August_LBCEQ_distance_vector,August_data))
summary(lm(August_unweighted_unifrac_distance_vector ~ August_LBCEQ_distance_vector,August_data))
summary(lm(August_bray_curtis_distance_vector ~ August_LBCEQ_distance_vector,August_data))


#########################################
library(ggpmisc)
#plot two vector 
###
August_N <- ggplot(data=August_data, aes(x=August_N_distance_vector, y=August_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in soil N level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 25),
    legend.title = element_text(color = "black", size = 25),
    axis.title.x = element_text(size=25, face="bold"),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=25, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(0.042,0.04),
           y = c(0.13,0.12), 
           label = c("p-value = 6.34e-08","R ^ 2 = 0.3592"),size=10) 
August_N

ggsave("August_N_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("August_N_weighted_unifrac_distance.png", height=8, width=10, device="png")

August_LBCEQ <- ggplot(data=data, aes(x=LBCEQ_distance_vector, y=unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in LBCEQ level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 25),
    legend.title = element_text(color = "black", size = 25),
    axis.title.x = element_text(size=25, face="bold"),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=25, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(160,160),
           y = c(0.13,0.12), 
           label = c("p-value = 1.817e-05","R ^ 2 = 0.2395"),size=10) 

August_LBCEQ 

ggsave("August_LBCEQ_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("August_LBCEQ_weighted_unifrac_distance.png", height=8, width=10, device="png")

###########################################################

## read in august weighted unifrac distance 
may_weighted_unifrac_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-weighted-unifrac-distance-matrix.tsv",sep="\t",header=TRUE)
may_weighted_unifrac_distance_vector <- calc_distance_vector(may_weighted_unifrac_distance)
may_weighted_unifrac_distance_vector

##
may_unweighted_unifrac_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-unweighted-unifrac-distance-matrix.tsv",sep="\t",header=TRUE)
may_unweighted_unifrac_distance_vector <- calc_distance_vector(may_unweighted_unifrac_distance)
may_unweighted_unifrac_distance_vector

may_bray_curtis_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-bray-curtis-distance-matrix.tsv",sep="\t",header=TRUE)
may_bray_curtis_distance_vector <- calc_distance_vector(may_bray_curtis_distance)
may_bray_curtis_distance_vector

length(may_bray_curtis_distance_vector)

## read in august N unifrac distance 
may_Zn_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/May-Zn-distance-matrix.tsv",sep="\t",header=TRUE)
##delete previous sample names column
rownames(may_Zn_distance) <- may_Zn_distance$X
may_Zn_distance$X <-NULL
may_Zn_distance <-may_Zn_distance[order(row.names(may_Zn_distance)),  order(colnames(may_Zn_distance))]
may_Zn_distance <- may_Zn_distance[-c(3),-c(3)]

ncol(may_Zn_distance)

#set up loop numbers 
row_loop=11
column_loop=10
##set up an emptry vector
may_Zn_distance_vector <- c(0)

######

###### loop trhough the lower triangle of the matrix and add elements to the vector 
while(row_loop <=11 & row_loop >1){
  for(i in column_loop:1){
    print(i)
    print(May_Zn_distance[row_loop,i])
    may_Zn_distance_vector <- append(may_Zn_distance_vector,may_Zn_distance[row_loop,i])
  }
  column_loop=column_loop-1
  row_loop = row_loop-1
}
### remove first zero 
may_Zn_distance_vector <- may_Zn_distance_vector[2:length(may_Zn_distance_vector)]
may_Zn_distance_vector
length(may_Zn_distance_vector)


############################

May_data <- data.frame(may_weighted_unifrac_distance_vector,may_unweighted_unifrac_distance_vector,may_bray_curtis_distance_vector,may_Zn_distance_vector)
summary(lm(may_weighted_unifrac_distance_vector ~ may_Zn_distance_vector, data=May_data))
summary(lm(may_unweighted_unifrac_distance_vector ~ may_Zn_distance_vector, data=May_data))
summary(lm(may_bray_curtis_distance_vector ~ may_Zn_distance_vector, data=May_data))



###################################################################################################

ggplot(data=May_data, aes(x=may_Zn_distance_vector, y=may_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in Zn level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 25),
    legend.title = element_text(color = "black", size = 25),
    axis.title.x = element_text(size=25, face="bold"),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=25, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(1.75,1.75),
           y = c(0.13,0.12), 
           label = c("p-value = 5.293e-06","R ^ 2 = 0.3134"),size=10) 

ggsave("Zn_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("Zn_weighted_unifrac_distance.png", height=8, width=10, device="png")


################

## read in august weighted unifrac distance 
June_weighted_unifrac_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-weighted-distance-matrix.tsv",sep="\t",header=TRUE)
June_weighted_unifrac_distance_vector <- calc_distance_vector(June_weighted_unifrac_distance) 

June_unweighted_unifrac_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-unweighted-distance-matrix.tsv",sep="\t",header=TRUE)
June_unweighted_unifrac_distance_vector <- calc_distance_vector(June_unweighted_unifrac_distance) 

June_bray_curtis_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-bray-curtis-distance-matrix.tsv",sep="\t",header=TRUE)
June_bray_curtis_distance_vector <- calc_distance_vector(June_bray_curtis_distance) 

## read in august N unifrac distance 
June_N_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/N-June-distance-matrix.tsv",sep="\t",header=TRUE)
June_N_distance_vector <- calc_distance_vector(June_N_distance) 

## read in august N unifrac distance 
June_LBCEQ_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-LBCEQ-distance-matrix.tsv",sep="\t",header=TRUE)
June_LBCEQ_distance_vector <- calc_distance_vector(June_LBCEQ_distance) 

#########################################################

June_data <- data.frame(June_N_distance_vector,June_LBCEQ_distance_vector,June_weighted_unifrac_distance_vector,June_unweighted_unifrac_distance_vector,June_bray_curtis_distance_vector)

summary(lm(June_weighted_unifrac_distance_vector ~ June_N_distance_vector,data=June_data))
summary(lm(June_unweighted_unifrac_distance_vector ~ June_N_distance_vector,data=June_data))
summary(lm(June_bray_curtis_distance_vector ~ June_N_distance_vector,June_data))

summary(lm(June_weighted_unifrac_distance_vector ~ June_LBCEQ_distance_vector,June_data))
summary(lm(June_unweighted_unifrac_distance_vector ~ June_LBCEQ_distance_vector,June_data))
summary(lm(June_bray_curtis_distance_vector ~ June_LBCEQ_distance_vector,June_data))



####################################################

data <- data.frame(June_N_distance_vector,June_unifrac_distance_vector,June_LBCEQ_distance_vector)

June_N <- ggplot(data=data, aes(x=June_N_distance_vector, y=June_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in soil N level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 25),
    legend.title = element_text(color = "black", size = 25),
    axis.title.x = element_text(size=25, face="bold"),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=25, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(0.05,0.05),
           y = c(0.13,0.12), 
           label = c("p-value = 0.001835","R ^ 2 = 0.1283"),size=10) 

June_N
ggsave("June_N_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("June_N_unifrac_distance.png", height=8, width=10, device="png")


June_LBCEQ <- ggplot(data=data, aes(x=June_LBCEQ_distance_vector, y=June_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  xlab("distance in LBCEQ level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 25),
    legend.title = element_text(color = "black", size = 25),
    axis.title.x = element_text(size=25, face="bold"),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=25, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(200,200),
           y = c(0.14,0.13), 
           label = c("p-value = 3.741e-05","R ^ 2 = 0.2229"),size=10) 

June_LBCEQ

ggsave("June_LBCEQ_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("June_LBCEQ_unifrac_distance.png", height=8, width=10, device="png")


#############################################
meta_soil_beta <- ggarrange(June_N,August_N,June_LBCEQ,August_LBCEQ,
          labels = c("June","August","June","August"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
          vjust = 2,
          hjust=0,
          ncol = 2, 
          nrow = 2)
meta_soil_beta 

ggsave("soil_physical_data_correlation_weighted_unifrac_distance.pdf",height=16, width=20, device="pdf")
ggsave("soil_physical_data_correlation_weighted_unifrac_distance.png",height=16, width=20, device="png")

aes(color="blue")
?annotate()
?label()

summary(lm(may_unifrac_distance_vector ~ May_Zn_distance_vector, data=data))
summary(lm(June_unifrac_distance_vector ~ June_N_distance_vector))
summary(lm(June_unifrac_distance_vector ~ June_LBCEQ_distance_vector))

citation("emmeans")
citation("car")


##########################################
#make supplement graph for 
sample_frequency <- read.table("/home/hl46161/Downloads/sample-frequency-detail(2).csv",sep=",")
colnames(sample_frequency) <- c("SampleID","Frequency")
?geom_histogram()
ggplot(data=sample_frequency,aes(Frequency)) + geom_histogram(fill='blue',col="red", alpha=0.2) +  
  ggtitle("                        Histogram for sample reads size distribution") +
  xlab('Number of reads') +
  ylab("Number of Sample") +
  theme(
    legend.text = element_text(color = "black", size = 18),
    legend.title = element_text(color = "black", size = 18),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"),
    axis.text.y = element_text(size=18, face="bold"),
    axis.text.x = element_text(size=18, face="bold"),
    plot.title = element_text(size=18, face="bold")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("reads_distribution.pdf", height=8, width=10, device="pdf")
ggsave("reads_distribution.png", height=8, width=10, device="png")


########################################### make rarefraction graph
install.packages ("rjson")
library("reshape2")
library("rjson")
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)

### read in the ASV count file 
observed_features_rarefration_data <- read.csv("/home/hl46161/publish_living_mulch/observed_features.csv")
### summarize the file based on sample ID
observed_features_rarefration_data_melt <- melt(observed_features_rarefration_data,id.vars = "sample.id")
### order the file based on sample ID 
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[order(observed_features_rarefration_data_melt$sample.id),]
n <- 10
##### summarize mean value of AVS count of every 10 iteration at one sequence septh
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt %>% group_by(mean = (row_number() -1) %/% n) %>%
  mutate(mean = mean(value))

#### sustitute depth.1_iter.1 to contain only sequence depth 
observed_features_rarefration_data_melt$variable <- gsub("depth.","",observed_features_rarefration_data_melt$variable)
observed_features_rarefration_data_melt$variable <- gsub("_iter.10","",observed_features_rarefration_data_melt$variable)
observed_features_rarefration_data_melt$variable <- gsub("_iter.[0-9]","",observed_features_rarefration_data_melt$variable)
##### reselect the column 
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[,c("sample.id","variable","mean")]
##### remove duplicate rows 
observed_features_rarefration_data_melt <- unique(observed_features_rarefration_data_melt)
##### rename the colulmn 
colnames(observed_features_rarefration_data_melt) <- c("SampleID","depth","ASV_count")
########
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[order(observed_features_rarefration_data_melt$depth),]
###### set the depth to bu numeric. this step is important. you can try not to do it to see what happen to the graph
observed_features_rarefration_data_melt$depth <- as.numeric(observed_features_rarefration_data_melt$depth)
ggplot(data=observed_features_rarefration_data_melt,aes(x =depth, y = ASV_count,group=SampleID,colour=SampleID)) + geom_line() + geom_point() +
  ylim(0, 5000) +
  scale_x_continuous("read depth",breaks=c(0,5000,10000,15000,20000,25000,30000),limits = c(0,32000)) +
  geom_text(x=25000, y=600, label="GA18-20S-600",color="brown") +
  theme(legend.position='none',axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text = element_text(size=15, face="bold")
  )

ggsave("rarefraction_curve_observed_ASV.pdf",height=10, width=8, device="pdf")
ggsave("rarefraction_curve_observed_ASV.png",height=10, width=8, device="png")


#############################################################################

### read in the shannon count file 
shannon_rarefration_data <- read.csv("/home/hl46161/publish_living_mulch/shannon.csv")
### summarize the file based on sample ID
shannon_rarefration_data_melt <- melt(shannon_rarefration_data,id.vars = "sample.id")
### order the file based on sample ID 
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[order(shannon_rarefration_data_melt$sample.id),]
n <- 10
##### summarize mean value of AVS count of every 10 iteration at one sequence septh
shannon_rarefration_data_melt <- shannon_rarefration_data_melt %>% group_by(mean = (row_number() -1) %/% n) %>%
  mutate(mean = mean(value))

#### sustitute depth.1_iter.1 to contain only sequence depth 
shannon_rarefration_data_melt$variable <- gsub("depth.","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt$variable <- gsub("_iter.10","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt$variable <- gsub("_iter.[0-9]","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[,c("sample.id","variable","mean")]
shannon_rarefration_data_melt <- unique(shannon_rarefration_data_melt)
##### rename the colulmn 
colnames(shannon_rarefration_data_melt) <- c("SampleID","depth","Shannon_index")
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[order(shannon_rarefration_data_melt$depth),]
###### set the depth to bu numeric. this step is important. you can try not to do it to see what happen to the graph
shannon_rarefration_data_melt$depth <- as.numeric(shannon_rarefration_data_melt$depth)
ggplot(data=shannon_rarefration_data_melt,aes(x =depth, y = Shannon_index,group=SampleID,color=SampleID)) + geom_line() + geom_point() + 
  ylim(0, 12) + 
  geom_text(x=25000, y=7.8, label="GA18-20S-600",color="brown") +
  scale_x_continuous("read depth",breaks=c(0,5000,10000,15000,20000,25000,30000),limits = c(0,32000)) +
  theme(legend.position='none',axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text = element_text(size=15, face="bold")
  )
ggsave("rarefraction_curve_shannon.pdf",height=10, width=8, device="pdf")
ggsave("rarefraction_curve_shannon.png",height=10, width=8, device="png")

###########################################
library(vegan)
#### the txt filed need to remove the # created from biom header before imported into R
otu_table <- read.table("/home/hl46161/publish_living_mulch/exported-feature-table-for-NMDS/living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000-rarefied.txt",sep="\t",header=TRUE,row.names = 1,check.names = FALSE)
otu_table <- t(otu_table)
otu_table <- as.data.frame(otu_table)
otu_table$SampleID <- rownames(otu_table)
otu_table$SampleID

metadata <- read.table("/home/hl46161/publish_living_mulch/living_mulch_data_metafile_massalin.tsv",sep="\t",header=TRUE)
merge_file <- merge(metadata,otu_table,by="SampleID")

August_samples <- merge_file[merge_file$Date_taken == "2018-08-31",]
June_samples <- merge_file[merge_file$Date_taken == "2018-06-28",]
May_samples <- merge_file[merge_file$Date_taken == "2018-05-21",]

August_abund <- August_samples[,c(20:ncol(August_samples))]
June_abund <- June_samples[,c(20:ncol(June_samples))]
May_abund <- May_samples[,c(20:ncol(May_samples))]

August_soil_physical_value <- August_samples[,c(5:18)]
June_soil_physical_value <- June_samples[,c(5:18)]
May_soil_physical_value <- May_samples[,c(5:18)]

#############
## run NMDS based on bray curtis distance 
May_nmds = metaMDS(May_abund, distance = "bray")
May_nmds
## fit  enviromental varibales into NMDS
en = envfit(May_nmds,May_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(May_nmds)
plot(en)

May_nmds$points

data.scores = as.data.frame(scores(May_nmds)) ## extract NMDS information 
data.scores$treatment = May_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

may_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + 
  scale_shape_manual(values=c(15, 16, 17,18)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 20, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 20, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  labs(colour = "Treatment",shape ="Treatment")



may_gg

ggsave("/home/hl46161/publish_living_mulch/May_NMDS.pdf", height=8, width=8, device="pdf")

###############################################################################################
## run NMDS based on bray curtis distance 
June_nmds = metaMDS(June_abund, distance = "bray")
June_nmds
## fit  enviromental varibales into NMDS
en = envfit(June_nmds,June_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(June_nmds)
plot(en)

June_nmds$points

data.scores = as.data.frame(scores(June_nmds)) ## extract NMDS information 
data.scores$treatment = June_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.045) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

June_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + 
  scale_shape_manual(values=c(15, 16, 17,18)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 20, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 20, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  labs(colour = "Treatment",shape ="Treatment")

June_gg

ggsave("/home/hl46161/publish_living_mulch/June_NMDS.pdf", height=8, width=8, device="pdf")

#,legend.position='none'

###########################################################################################
## run NMDS based on bray curtis distance 
August_nmds = metaMDS(August_abund, distance = "bray")
August_nmds
## fit  enviromental varibales into NMDS
en = envfit(August_nmds,August_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(August_nmds)
plot(en)

August_nmds$points

data.scores = as.data.frame(scores(August_nmds)) ## extract NMDS information 
data.scores$treatment = August_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

August_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + 
  scale_shape_manual(values=c(15, 16, 17,18)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 20, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 20, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  labs(colour = "Treatment",shape ="Treatment")

August_gg

ggsave("/home/hl46161/publish_living_mulch/August_NMDS.pdf", height=8, width=8, device="pdf")

########################################################
library(gridExtra)
library(ggpubr)

mega_NMDS = ggarrange(may_gg,June_gg,August_gg,
                      labels = c("May","June", "August"),
                      font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                      vjust = 5,
                      hjust=0,
                      ncol = 3, 
                      nrow = 1,
                      common.legend = TRUE)

?ggarrange()
mega_NMDS

ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.pdf", height=12, width=25, device="pdf")
ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.png", height=12, width=25, device="png")






