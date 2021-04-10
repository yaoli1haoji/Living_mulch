library("phyloseq")
library("ggplot2")
#install.packages("ggplot")
library("vegan")
library("DESeq2")
#library("pathview")
#install.packages("RCurl")
#install.packages("igraph")
#library(igraph)
library("RCy3")
library("dplyr")
#install.packages("dplyr")
#install.packages("stringgaussnet")
library(microbiome)
library(ggnet)
install.packages("remotes")
remotes::install_github("steinbaugh/DESeqAnalysis")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")

ainstalled.packages()[, c("Package", "LibPath")]
BiocManager::install("RCy3")
aaif(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
Sys.setenv(R_INSTALL_STAGED = FALSE)
library(RCy3)
library(igraph)

source('../utility/cytoscape_util.R')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("viridis")
library(viridis)

devtools::install_github("briatte/ggnet")
library(ggnet)


source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
browseVignettes("DESeq2")

#setwd("/home/hl46161/new_investigation_living_mulch/")

#import otu table and convert to matrix. This is the otu table after mitochondira and chloroplast filtration, no blank
otu_table	=	read.csv("/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.txt",sep="\t",row.names =1,check.names = FALSE)

#import otu table and convert to matrix 
otu_table	=	as.matrix(otu_table)

#taxonomy come from qiime2 artifact and needs to seperate in domain, phylum,order,......
taxonomy = read.csv("/home/hl46161/publish_living_mulch/taxonomy.tsv",sep="\t",row.names =1)
#taxonomy = read.csv("taxonomy.tsv",sep="\t",row.names =1)
taxonomy	=	as.matrix(taxonomy)
#import phylogenetic tree table and convert to matrix
phy_tree	=	read_tree("/home/hl46161/publish_living_mulch/tree.nwk")
#import metadata table and convert to matrix 
metadata	=	read.table("/home/hl46161/publish_living_mulch/living_mulch_data_metafile_massalin.tsv",sep ="\t",header = 1,row.names=1)
#metadata <- metadata[-c(4),]

#convert the data into phyloseq format 
OTU	=	otu_table(otu_table,taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META = sample_data(metadata)
#make sure taxonomy and otu table match. The differences come from filtered mitochondira and cloroplast 
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

sample_names(OTU)
sample_names(META)

##build phyloseq object 
ps	=	phyloseq(OTU,TAX,phy_tree,META)
## create rarefaction plot to see the read depth for each samples 
rarecurve(t(otu_table(ps)), step=50, cex=0.5)
## save the rarefraction graph if necessary 
ggsave("rarefraction")


## create rarefied phyloseq object  size = minimium sample size * 0.9 
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)
## seperate data set into different time points 
ps.mid = subset_samples(ps, Date_taken == "2018-06-28")
ps.early = subset_samples(ps, Date_taken == "2018-05-21")
ps.late = subset_samples(ps, Date_taken == "2018-08-31")

ps.lm = subset_samples(ps, treatment == "LivingMulch")

ps.late.rarefied = subset_samples(ps.rarefied, Date_taken == "2018-08-31")
sample_sums(ps)
sample_sums(ps.rarefied)
plot_bar(ps,fill = "phylum")

otu_table(ps)
#######################################
##alpha diversity
plot_richness(ps.rarefied, x="treatment", color="Date_taken", measures=c("Observed"))

plot_richness(ps.rarefied, x="Date_taken", color="treatment", measures=c("Observed","Shannon")) + geom_boxplot()

plot_richness(ps.rarefied, x="treatment", measures=c("Observed", "Shannon"),color = "Date_taken") + geom_boxplot()

#beta diversity test
rich = estimate_richness(ps.rarefied)
rich
pairwise.wilcox.test(rich$Observed, sample_data(ps.rarefied)$treatment)
pairwise.wilcox.test(rich$Shannon, sample_data(ps.rarefied)$treatment)
plot_richness(ps.rarefied, x="Date_taken",color="treatment", measures=c("Observed", "Shannon")) + geom_boxplot()

ps.late.rarefied <- subset_samples(ps.rarefied, Date_taken = "2018-08-31")
rich2 = estimate_richness(ps.late.rarefied)
pairwise.wilcox.test(rich2$Observed, sample_data(ps.late.rarefied)$treatment)
pairwise.wilcox.test(rich2$Shannon, sample_data(ps.late.rarefied)$treatment)

rich3 = estimate_richness(ps)
pairwise.wilcox.test(rich3$Observed, sample_data(ps)$treatment)
pairwise.wilcox.test(rich3$Shannon, sample_data(ps)$treatment)

##PCOA graph construction 
wunifrac_dist = phyloseq::distance(ps.late.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps.late.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="treatment") + theme(aspect.ratio=1)
adonis(wunifrac_dist ~ sample_data(ps.rarefied)$treatment)


#######################################################################################

#differential abundance test 
#import phyloseq oject into deseq2
ds = phyloseq_to_deseq2(ps, ~ Treatment)
ds = DESeq(ds)
#set the detection threshold 
alpha = 0.05

##first calculate overall differential abundant OTU between No Cover and living mulch 
res = results(ds, contrast=c("Treatment","LivingMulch","NoCover"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
#merge the otu with taxonomy
lm_vs_nc_ps = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
#filter based on critieria
lm_vs_nc_filtered_1_ps <- subset(lm_vs_nc_ps, abs(log2FoldChange) > 0)
lm_vs_nc_filtered_1_ps
lm_vs_nc_filtered_1_ps <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_vs_nc_filtered_1_ps

#output the results of lm vs nc 
write.csv(lm_vs_nc_filtered_1_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/lm_VS_nocover_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and living mulch 
res2 = results(ds, contrast=c("Treatment","LivingMulch","CrimsonClover"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2


##first calculate overall differential abundant OTU between CerealRye and living mulch 
res4 = results(ds, contrast=c("Treatment","LivingMulch","CerealRye"), alpha=alpha)
res4 = res4[order(res4$padj, na.last=NA), ]
res_sig4 = res4[(res4$padj < alpha), ]
res_sig4
lm_vs_cr_4_ps = cbind(as(res_sig4, "data.frame"), as(tax_table(ps)[rownames(res_sig4), ], "matrix"))
lm_vs_cr_filtered_4_ps <- subset(lm_vs_cr_4_ps, baseMean > 0)
lm_vs_cr_filtered_4_ps
lm_vs_cr_filtered_4_ps <- subset(lm_vs_cr_filtered_4_ps,abs(log2FoldChange) > 0)
lm_vs_cr_filtered_4_ps

write.csv(lm_vs_cr_filtered_4_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/lm_VS_cr_differential_OTU.csv")

##first calculate overall differential abundant OTU between CerealRye and Nocover 

res5 = results(ds, contrast=c("Treatment", "CerealRye","NoCover"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig5 = res5[(res5$padj < alpha), ]
res_sig5
cr_vs_nc_5_ps = cbind(as(res_sig5, "data.frame"), as(tax_table(ps)[rownames(res_sig5), ], "matrix"))
cr_vs_nc__filtered_5_ps <- subset(cr_vs_nc_5_ps, baseMean > 0)
cr_vs_nc__filtered_5_ps
cr_vs_nc__filtered_5_ps <- subset(cr_vs_nc__filtered_5_ps,abs(log2FoldChange) > 0)
cr_vs_nc__filtered_5_ps

write.csv(cr_vs_nc__filtered_5_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/cr_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res6 = results(ds, contrast=c("Treatment", "CrimsonClover","NoCover"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig6 = res6[(res6$padj < alpha), ]
res_sig6
cc_vs_nc_6 = cbind(as(res_sig6, "data.frame"), as(tax_table(ps)[rownames(res_sig6), ], "matrix"))
cc_vs_nc__filtered_6_ps <- subset(cc_vs_nc_6, baseMean > 0)
cc_vs_nc__filtered_6_ps
cc_vs_nc__filtered_6_ps <- subset(cc_vs_nc__filtered_6_ps,abs(log2FoldChange) > 0)
cc_vs_nc__filtered_6_ps

write.csv(cc_vs_nc__filtered_6_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/cc_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res7 = results(ds, contrast=c("Treatment", "CrimsonClover","CerealRye"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig7 = res7[(res7$padj < alpha), ]
res_sig7
cc_vs_cr_7_ps = cbind(as(res_sig7, "data.frame"), as(tax_table(ps)[rownames(res_sig7), ], "matrix"))
cc_vs_cr__filtered_7_ps <- subset(cc_vs_cr_7_ps, baseMean > 0)
cc_vs_cr__filtered_7_ps
cc_vs_cr__filtered_7_ps <- subset(cc_vs_cr__filtered_7_ps,abs(log2FoldChange) > 0)
cc_vs_cr__filtered_7_ps 

write.csv(cc_vs_cr__filtered_7_ps,"/home/hl46161/living_mulch_result_summry/CC_VS_CR_differential_OTU.csv")


#generate a ven diagram
gg <- gplots::venn(list(LivingMulch=rownames(lm_vs_nc_ps),CerealRye = rownames(cc_vs_nc_6),CrimsonColver = rownames(cr_vs_nc_5_ps)))
ggsave(gg,"Venn_diagram", height=6, width=6, device="pdf")

####################################################################

ps.lm
ds2 = phyloseq_to_deseq2(ps.lm, ~Date_taken)
ds2 = DESeq(ds2)
alpha = 0.05

res = results(ds2, contrast=c("Date_taken","2018-05-21","2018-08-31"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
lm_may_vs_Au = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
lm_may_vs_Au  <- subset(lm_may_vs_Au, abs(log2FoldChange) > 0)
lm_may_vs_Au 
lm_may_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_may_vs_Au 


res2 = results(ds2, contrast=c("Date_taken","2018-06-28","2018-08-31"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_june_vs_Au = cbind(as(res_sig2, "data.frame"), as(tax_table(ps)[rownames(res_sig2), ], "matrix"))
lm_june_vs_Au  <- subset(lm_june_vs_Au, abs(log2FoldChange) > 0)
lm_june_vs_Au 
lm_june_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_june_vs_Au 



#####################################################################


ds = phyloseq_to_deseq2(ps, ~ Date_taken)
ds = DESeq(ds)
otu_table(ps)
alpha = 0.05

res = results(ds2, contrast=c("Date_taken","2018-05-21","2018-08-31"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
may_vs_Au = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
may_vs_Au  <- subset(may_vs_Au, abs(log2FoldChange) > 0)
may_vs_Au 
may_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
may_vs_Au 


res2 = results(ds2, contrast=c("Date_taken","2018-06-28","2018-08-31"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_june_vs_Au = cbind(as(res_sig2, "data.frame"), as(tax_table(ps)[rownames(res_sig2), ], "matrix"))
lm_june_vs_Au  <- subset(lm_june_vs_Au, abs(log2FoldChange) > 0)
lm_june_vs_Au 
lm_june_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_june_vs_Au 




#####################################################calculate the differential abundant OTU based on taxnomic level 
?merge_taxa()

OrderFiltered <- subset_taxa(ps, Order != "NA")
ps.order = tax_glom(OrderFiltered , "Order")
otu_table(OrderFiltered)
otu_table(ps.order)
order_otu_table = as(otu_table(ps.order), "matrix")
order_otu_table  = as.data.frame(order_otu_table)
write.csv(order_otu_table,"/home/hl46161/publish_replicate_new_living_mulch/order_otu_table.csv")

FamilyFiltered <- subset_taxa(ps, Family != "NA")
ps.family = tax_glom(FamilyFiltered, "Family")
otu_table(FamilyFiltered)
otu_table(ps.family)
family_otu_table = as(otu_table(ps.family), "matrix")
family_otu_table  = as.data.frame(family_otu_table)
#write.table(family_otu_table,"/home/hl46161/publish_replicate_new_living_mulch/family_otu_table.tsv",sep="\t")

#differential abundance test 

ds = phyloseq_to_deseq2(ps.family, ~ Treatment)
ds = DESeq(ds)
alpha = 0.05

##first calculate overall differential abundant OTU between No Cover and living mulch 
res = results(ds, contrast=c("Treatment","LivingMulch","NoCover"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
lm_vs_nc_ps.family = cbind(as(res_sig, "data.frame"), as(tax_table(ps.family)[rownames(res_sig), ], "matrix"))
lm_vs_nc_filtered_1_ps.family <- subset(lm_vs_nc_ps.family, abs(log2FoldChange) > 0)
lm_vs_nc_filtered_1_ps.family
lm_vs_nc_filtered_1_ps.family <- subset(lm_vs_nc_filtered_1_ps.family, baseMean > 0)


write.csv(lm_vs_nc_filtered_1_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_nc_differential_OTU.csv")


##calculate overall differential abundant OTU between CrimsonClover and living mulch 
res2 = results(ds, contrast=c("Treatment","LivingMulch","CrimsonClover"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_vs_cc_2_ps.family = cbind(as(res_sig2, "data.frame"), as(tax_table(ps.family)[rownames(res_sig2), ], "matrix"))
lm_vs_cc_2_ps.family <- subset(lm_vs_cc_2_ps.family, baseMean > 0)
lm_vs_cc_2_ps.family
lm_vs_cc_2_ps.family <- subset(lm_vs_cc_2_ps.family,abs(log2FoldChange) > 0)
lm_vs_cc_2_ps.family


write.csv(lm_vs_cc_2_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CerealRye and living mulch 
res4 = results(ds, contrast=c("Treatment","LivingMulch","CerealRye"), alpha=alpha)
res4 = res4[order(res4$padj, na.last=NA), ]
res_sig4 = res4[(res4$padj < alpha), ]
res_sig4
lm_vs_cr_4_ps.family = cbind(as(res_sig4, "data.frame"), as(tax_table(ps.family)[rownames(res_sig4), ], "matrix"))
lm_vs_cr_filtered_4_ps.family <- subset(lm_vs_cr_4_ps.family, baseMean > 0)
lm_vs_cr_filtered_4_ps.family
lm_vs_cr_filtered_4_ps.family <- subset(lm_vs_cr_filtered_4_ps.family,abs(log2FoldChange) > 0)
lm_vs_cr_filtered_4_ps.family

write.csv(lm_vs_cr_filtered_4_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cr_differential_OTU.csv")


##first calculate overall differential abundant OTU between CerealRye and Nocover

res5 = results(ds, contrast=c("Treatment", "CerealRye","NoCover"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig5 = res5[(res5$padj < alpha), ]
res_sig5
cr_vs_nc_5_ps.family = cbind(as(res_sig5, "data.frame"), as(tax_table(ps.family)[rownames(res_sig5), ], "matrix"))
cr_vs_nc__filtered_5_ps.family <- subset(cr_vs_nc_5_ps.family, baseMean > 0)
cr_vs_nc__filtered_5_ps.family
cr_vs_nc__filtered_5_ps.family <- subset(cr_vs_nc__filtered_5_ps.family,abs(log2FoldChange) > 0)
cr_vs_nc__filtered_5_ps.family

write.csv(cr_vs_nc__filtered_5_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cr_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res6 = results(ds, contrast=c("Treatment", "CrimsonClover","NoCover"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig6 = res6[(res6$padj < alpha), ]
res_sig6
cc_vs_nc_6 = cbind(as(res_sig6, "data.frame"), as(tax_table(ps.family)[rownames(res_sig6), ], "matrix"))
cc_vs_nc__filtered_6_ps.family <- subset(cc_vs_nc_6, baseMean > 0)
cc_vs_nc__filtered_6_ps.family
cc_vs_nc__filtered_6_ps.family <- subset(cc_vs_nc__filtered_6_ps.family,abs(log2FoldChange) > 0)
cc_vs_nc__filtered_6_ps.family 

write.csv(cc_vs_nc__filtered_6_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res7 = results(ds, contrast=c("Treatment", "CrimsonClover","CerealRye"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig7 = res7[(res7$padj < alpha), ]
res_sig7
cc_vs_cr_7_ps.family = cbind(as(res_sig7, "data.frame"), as(tax_table(ps.family)[rownames(res_sig7), ], "matrix"))
cc_vs_cr__filtered_7_ps.family <- subset(cc_vs_cr_7_ps.family, baseMean > 0)
cc_vs_cr__filtered_7_ps.family
cc_vs_cr__filtered_7_ps.family <- subset(cc_vs_cr__filtered_7_ps.family,abs(log2FoldChange) > 0)
cc_vs_cr__filtered_7_ps.family 

write.csv(cc_vs_cr__filtered_7_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_cr_differential_OTU.csv")

#######################

#summarize differential abudant otu. First try to find out which dfOtu is shared between three treatment 
lm_vs_nc_filtered_1_ps.family$OTU_ID <- rownames(lm_vs_nc_filtered_1_ps.family)
cr_vs_nc__filtered_5_ps.family$OTU_ID <- rownames(cr_vs_nc__filtered_5_ps.family)
cc_vs_nc__filtered_6_ps.family$OTU_ID <- rownames(cc_vs_nc__filtered_6_ps.family)

lm_cr_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family,cc_vs_nc__filtered_6_ps.family,by = "OTU_ID")
lm_cr_cc_otu <- lm_cr_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,baseMean.x,padj.x)
lm_cr_cc_otu <- dplyr::inner_join(lm_cr_cc_otu, cr_vs_nc__filtered_5_ps.family,by = "OTU_ID")

lm_cr_cc_otu <- lm_cr_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)

lm_cr_cc_otu_100 <- subset(lm_cr_cc_otu,baseMean.x>=100)
lm_cr_cc_otu_100 <- subset(lm_cr_cc_otu,log2FoldChange.x>0)

#lm_cr_cc_otu <- lm_cr_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,log2FoldChange,baseMean.x,padj.x)
colnames(lm_cr_cc_otu) <-c("Kingdom","Phylum","Class","Order","Family","OTU_ID","log2FoldChange.lm","log2FoldChange.cc","log2FoldChange.cr","baseMean","P_adjusted")


lm_cr_cc_otu$log2FoldChange.lm <- round(lm_cr_cc_otu$log2FoldChange.lm, digits = 2)
lm_cr_cc_otu$log2FoldChange.cc <- round(lm_cr_cc_otu$log2FoldChange.cc, digits = 2)
lm_cr_cc_otu$log2FoldChange.cr <- round(lm_cr_cc_otu$log2FoldChange.cr, digits = 2)
#lm_cr_cc_otu$P_adjusted <- round(lm_cr_cc_otu$P_adjusted, digits = 2)

#lm_cr_cc_otu %>%
#mutate(Taxa = if_else( grepl("[[:digit:]]", lm_cr_cc_otu$Order),Order,Family))
#grepl("[[:digit:]]", lm_cr_cc_otu$Family)
#gsub("[[:digit:]]",lm_cr_cc_otu$Order,lm_cr_cc_otu$Family, fixed = TRUE)
?gsub()

lm_cr_cc_otu$Taxa <- as.character(lm_cr_cc_otu$Family)
lm_cr_cc_otu$Taxa[c(2,3,18,19,22,25,26,27,28,29,30,31,32,34,36,37,38,43)] <- c("c__Ktedonobacteria","o__Ktedonobacterales","p__Chloroflexi","o__Tepidisphaerales",
                                                                                " o__Rokubacteriales"," c__Planctomycetes"," c__Ktedonobacteria"," p__WPS-2","c__Phycisphaerae	",
                                                                                "c__Ktedonobacteria","p__Acidobacteriota"," c__Gammaproteobacteria","o__Acidobacteriales",
                                                                                "o__Gaiellales","c__Vicinamibacteria","p__Acidobacteriota","o__Elsterales",
                                                                               "c__Anaerolineae")



lm_cr_cc_otu_table <- lm_cr_cc_otu %>%
  mutate(Abudance =
           case_when(baseMean <= 100 ~ "*",
                     baseMean > 1000 ~ "***",
                     baseMean >100 ~ "**",
                     ))

lm_cr_cc_otu_table  <- lm_cr_cc_otu_table  %>% select(log2FoldChange.lm:log2FoldChange.cr,Taxa,Abudance)
 


#output the shared OTU between three cover crop treatment
write.csv(lm_cr_cc_otu,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_cc_cr_VS_nc_differential_OTU.csv",row.names = FALSE)

write.csv(lm_cr_cc_otu_table,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_table_lm_cr_cc_vs_nc.csv",row.names = FALSE)
  
  
lm_cr_cc_otu$CerealRye <- 1
lm_cr_cc_otu$LivingMulch <- 1
lm_cr_cc_otu$CrimsonClover <- 1

lm_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu <- dplyr::left_join(lm_cc_otu,cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu[rowSums(is.na(lm_cc_otu)) > 4,]
lm_cc_otu <- lm_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu$CerealRye <- 0
lm_cc_otu$LivingMulch <- 1
lm_cc_otu$CrimsonClover <- 1



lm_cr_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu <- dplyr::left_join(lm_cr_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu[rowSums(is.na(lm_cr_otu)) > 4,]
lm_cr_otu <- lm_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu$CerealRye <- 1
lm_cr_otu$LivingMulch <- 1
lm_cr_otu$CrimsonClover <- 0



cc_cr_otu <- dplyr::inner_join(cc_vs_nc__filtered_6_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu <- dplyr::left_join(cc_cr_otu,lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu[rowSums(is.na(cc_cr_otu)) > 4,]
#cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Damily.x,OTU_ID)
cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu$CerealRye <- 1
cc_cr_otu$LivingMulch <- 0
cc_cr_otu$CrimsonClover <- 1




lm_nc_otu <- dplyr::left_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_nc_otu <- dplyr::left_join(lm_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#lm_nc_otu <- lm_nc_otu [,c(7:ncol(lm_nc_otu))]
lm_nc_otu$CerealRye <- 0
lm_nc_otu$LivingMulch <- 1
lm_nc_otu$CrimsonClover <- 0

cc_nc_otu <- dplyr::left_join(cc_vs_nc__filtered_6_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_nc_otu <- dplyr::left_join(cc_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cc_nc_otu <- cc_nc_otu [,c(7:ncol(cc_nc_otu))]
cc_nc_otu$CerealRye <- 0
cc_nc_otu$LivingMulch <- 0
cc_nc_otu$CrimsonClover <- 1

cr_nc_otu <- dplyr::left_join(cr_vs_nc__filtered_5_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cr_nc_otu <- dplyr::left_join(cr_nc_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cr_nc_otu <- cr_nc_otu [,c(7:ncol(cr_nc_otu))]
cr_nc_otu$CerealRye <- 1
cr_nc_otu$LivingMulch <- 0
cr_nc_otu$CrimsonClover <- 0


otu_summary_matrix <- dplyr::bind_rows(lm_nc_otu,cc_nc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cr_nc_otu)

otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cc_cr_otu)
otu_summary_matrix$Confidence.x <- NULL

upset_data <- otu_summary_matrix %>% select(CerealRye:CrimsonClover)

upset(upset_data,
      point.size = 3.2, line.size = 2,
      mainbar.y.label = "Distrubution of differential family", sets.x.label = "Total # of differential family",
      text.scale =c(2, 2, 2, 2, 2, 2)
)

ggsave("differential_otu_upset_plot.pdf",width = 8,height = 8,device = "pdf")

grid.text(
  "@littlemissdata",
  x = 0.90,
  y = 0.05,
  gp = gpar(
    fontsize = 10,
    fontface = 10
  )
)


lm_cr_cc_otu %>% select(ingdom.x:Family.x,OTU_ID)
otu_summary_ploting_matrix <- otu_summary_matrix %>% select(OTU_ID:CrimsonClover)
rownames(otu_summary_ploting_matrix) <- otu_summary_ploting_matrix$OTU_ID
otu_summary_ploting_matrix$OTU_ID <- NULL
otu_summary_ploting_matrix$log2FoldChange.x <- NULL
otu_summary_ploting_matrix$baseMean.x <- NULL
#otu_summary_matrix <- otu_summary_matrix[,c(1:11)]

net = network(otu_summary_ploting_matrix, directed = FALSE)
taxa = otu_summary_matrix$Family.x
taxa = as.character(taxa)
taxa = as.factor(taxa)
levels(taxa)
net %v% "family" = as.character(taxa)
taxa
levels(taxa) 
nb.cols = 49
y = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nb.cols)
y
names(y)
names(y) = levels(taxa)
net %v% "expression" = ifelse(otu_summary_matrix$log2FoldChange.x >0,"increased abudances","decrease abudances")
gg <- ggnet2(net,color = "family",palette = y,shape = "expression",label = c("LivingMulch", "CrimsonClover","CerealRye"),legend.size = 6, legend.position = "bottom",node.size = 4,label.size = 3)  +
  theme(panel.background = element_rect(color = "grey"))
gg
ggsave(plot=gg,"differential_otu",width = 10,height = 10,device = "pdf")

?ggsave()

lm_nc_otu <- subset(lm_vs_nc_filtered_1_ps.family,lm_vs_nc_filtered_1_ps.family$OTU_ID == lm_cc_otu$OTU_ID)

lm_vs_nc_filtered_1_ps.family$OTU_ID == lm_cc_otu$OTU_ID


lm_vs_nc_filtered_1_ps.family[lm_vs_nc_filtered_1_ps.family$OTU %in% lm_nc_otu$OTU_ID,]


##################make heatmap for family level taxonomy 
# Transform count data using the variance stablilizing transform
ds
ds$
  otu_table(ds)
deseq2VST <- rlog(ds)
deseq2VST
deseq2ResDF <- as.data.frame(significant_otu_list)
# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.05, "Significant", NA)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$OTU <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 2,])
#deseq2ResDF[deseq2ResDF <= .05 & abs(deseq2ResDF$log2FoldChange) > 1,]
#sigGenes <- rownames(deseq2ResDF)
sigGenes
deseq2VST <- deseq2VST[deseq2VST$OTU %in% sigGenes,]


# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
#colnames(deseq2VST) <- c("Family","samples","value")
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("OTU"))
#colnames((deseq2VST_wide)) %in%
#metadata_2 <-  subset(metadata,treatment=="NoCover"|treatment=="LivingMulch")


head(deseq2VST_wide)
head(deseq2VST_long)
metadata_2$SampleID = rownames(metadata_2)
colnames(deseq2VST_long)[2] <- c("SampleID")
deseq2VST_long = dplyr::right_join(deseq2VST_long,metadata,by.x="variable",by.y="SampleID")
head(deseq2VST_long)
deseq2VST_long = deseq2VST_long[,c(1,2,3,5)]
head(deseq2VST_long)
# Now overwrite our original data frame with the long format
#deseq2VST <- melt(deseq2VST, id.vars=c("OTU"))

# Make a heatmap
heatmap <- ggplot(deseq2VST_long, aes(x=treatment, y=OTU, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt")  + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


























####################################################################
## set the aplha threslhod to be 0.01
ds2 = phyloseq_to_deseq2(ps.late, ~ treatment)
ds2 = DESeq(ds2)
alpha = 0.01
res5 = results(ds2, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig_late = res5[(res5$padj < alpha), ]
res_sig_late
res_sig_late = cbind(as(res_sig_late, "data.frame"), as(tax_table(ps)[rownames(res_sig_late), ], "matrix"))
res_sig_filtered_late <- subset(res_sig_late,abs(log2FoldChange) > 1)
res_sig_filtered_late <- subset(res_sig_filtered_late, baseMean > 10)
res_sig_filtered_late
write.table(res_sig_filtered_late, "/home/hl46161/new_living_mulch/phyloseq_data/differential_OTU_filtered_late", sep="\t")
sample_names(res_sig_filtered_late)
plot_bar(res_sig_filtered_late,facet_grid=~SampleType,taxa_are_rows =TRUE)
?read.table()

?plot_bar
ds3 = phyloseq_to_deseq2(ps.mid.rarefied, ~ treatment)
ds3 = DESeq(ds3)
alpha = 0.01
res6 = results(ds3, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig_mid = res6[(res6$padj < alpha), ]
res_sig_mid
res_sig_mid = cbind(as(res_sig_mid, "data.frame"), as(tax_table(ps)[rownames(res_sig_mid), ], "matrix"))
res_sig_filtered_mid <- subset(res_sig_mid,abs(log2FoldChange) > 1)
res_sig_filtered_mid <- subset(res_sig_filtered_mid, baseMean > 10)
res_sig_filtered_mid

ds4 = phyloseq_to_deseq2(ps.early.rarefied, ~ treatment)
ds4 = DESeq(ds4)
alpha = 0.01
res7 = results(ds4, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig_early = res7[(res7$padj < alpha), ]
res_sig_early
res_sig_early = cbind(as(res_sig_early, "data.frame"), as(tax_table(ps)[rownames(res_sig_early), ], "matrix"))
res_sig_filtered_early <- subset(res_sig_early,abs(log2FoldChange) > 1)
res_sig_filtered_early <- subset(res_sig_filtered_early, baseMean > 10)
res_sig_filtered_early


lm_relation_shannon <- lm(ds$shannon ~ds$treatment*ds$Date_taken,data = merge_file)

### increase the cutoff to see what happen 

ds2 = phyloseq_to_deseq2(ps.late, ~ treatment)
ds2 = DESeq(ds2)
alpha = 0.05
res5 = results(ds2, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig_late = res5[(res5$padj < alpha), ]
res_sig_late
#res_sig_late = cbind(as(res_sig_late, "data.frame"), as(tax_table(ps)[rownames(res_sig_late), ], "matrix"))
res_sig_filtered_late_0.05 <- subset(res_sig_late,abs(log2FoldChange) > 1)
res_sig_filtered_late_0.05 <- subset(res_sig_filtered_late_0.05, baseMean > 10)
res_sig_filtered_late_0.05

    plotPCA()

?plot_bar
ds3 = phyloseq_to_deseq2(ps.mid, ~ treatment)
ds3 = DESeq(ds3)
alpha = 0.05
res6 = results(ds3, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig_mid_0.05 = res6[(res6$padj < alpha), ]
res_sig_mid_0.05
res_sig_mid_0.05 = cbind(as(res_sig_mid_0.05, "data.frame"), as(tax_table(ps)[rownames(res_sig_mid_0.05), ], "matrix"))
res_sig_filtered_mid_0.05 <- subset(res_sig_mid_0.05,abs(log2FoldChange) > 1)
res_sig_filtered_mid_0.05 <- subset(res_sig_filtered_mid_0.05, baseMean > 10)
res_sig_filtered_mid_0.05

ds4 = phyloseq_to_deseq2(ps.early, ~ treatment)
ds4 = DESeq(ds4)
alpha = 0.05
res7 = results(ds4, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig_early_0.05 = res7[(res7$padj < alpha), ]
res_sig_early_0.05
res_sig_early_0.05 = cbind(as(res_sig_early_0.05, "data.frame"), as(tax_table(ps)[rownames(res_sig_early_0.05), ], "matrix"))
res_sig_filtered_early_0.05 <- subset(res_sig_early_0.05,abs(log2FoldChange) > 1)
res_sig_filtered_early_0.05 <- subset(res_sig_filtered_early_0.05, baseMean > 10)
res_sig_filtered_early_0.05


consistent_otu_2 <- merge(res_sig_filtered_early_0.05,res_sig_filtered_late_0.05,by ="row.names")
class(otuTable.late)
##otuTable.late <- subset(otuTable.late , conserved_otu %in% colnames(otuTable.late))
otuTable.late <- otuTable.late[,(colnames(otuTable) %in% conserved_otu)]
otuTable.late <- otuTable.late[c(1,14:)]
ggplot(otuTable.late,aes(x=, y=colnames))

conserved_otu <- consistent_otu_2$Row.names
conserved_otu <- c(conserved_otu)
conserved_otu

??phyloseq

## keeping trac of diferential abundance OTU in early data set 
early_OTU_ID <- data.frame(OTUID = row.names(res_sig_filtered_early_0.05))
mid_OTU_ID <- data.frame(OTUID = row.names(res_sig_filtered_mid_0.05))
late_OTU_ID <- data.frame(OTUID = row.names(res_sig_filtered_late_0.05))
differntial_OTU_ID <- rbind.data.frame(early_OTU_ID,mid_OTU_ID)
differntial_OTU_ID <-  rbind.data.frame(differntial_OTU_ID,late_OTU_ID)
differntial_OTU_ID <- distinct(differntial_OTU_ID)


#increase the cut off to 1 so that can track the change of OTU in three time point 
ds2 = phyloseq_to_deseq2(ps.late, ~ treatment)
ds2 = DESeq(ds2)
alpha = 0.99999
res5 = results(ds2, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig_late = res5[(res5$padj < alpha), ]
res_sig_late
res_sig_late = cbind(as(res_sig_late, "data.frame"), as(tax_table(ps)[rownames(res_sig_late), ], "matrix"))
res_sig_filtered_late_1.00 <- subset(res_sig_late,abs(log2FoldChange) > 0)
##res_sig_filtered_late_1.00 <- subset(res_sig_filtered_late_1.00, baseMean > 10)
res_sig_filtered_late_1.00
## filter out extra OTU to keep track of siginificant OTU
res_sig_filtered_late_1.00 <- res_sig_filtered_late_1.00[(rownames(res_sig_filtered_late_1.00) %in% differntial_OTU_ID$OTUID),]
res_sig_filtered_late_1.00$Date = "8-31-2018" 
write.table(res_sig_filtered_late, "/home/hl46161/new_living_mulch/phyloseq_data/differential_OTU_filtered_late", sep="\t")
?plot_bar


ds3 = phyloseq_to_deseq2(ps.mid, ~ treatment)
ds3 = DESeq(ds3)
alpha = 0.9999
res6 = results(ds3, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig_mid_1.00 = res6[(res6$padj < alpha), ]
res_sig_mid_1.00
res_sig_mid_1.00 = cbind(as(res_sig_mid_1.00, "data.frame"), as(tax_table(ps)[rownames(res_sig_mid_1.00), ], "matrix"))
res_sig_filtered_mid_1.00 <- subset(res_sig_mid_1.00,abs(log2FoldChange) > 0)
##res_sig_filtered_mid_1.00 <- subset(res_sig_filtered_mid_1.00, baseMean > 10)
res_sig_filtered_mid_1.00
res_sig_filtered_mid_1.00 <- res_sig_filtered_mid_1.00[(rownames(res_sig_filtered_mid_1.00) %in% differntial_OTU_ID$OTUID),]
res_sig_filtered_mid_1.00$Date = "6-21-2018" 

ds4 = phyloseq_to_deseq2(ps.early, ~ treatment)
ds4 = DESeq(ds4)
alpha = 0.99999
res7 = results(ds4, contrast=c("treatment", "NoCover","LivingMulch"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig_early_1.00 = res7[(res7$padj < alpha), ]
res_sig_early_1.00
res_sig_early_1.00 = cbind(as(res_sig_early_1.00, "data.frame"), as(tax_table(ps)[rownames(res_sig_early_1.00), ], "matrix"))
res_sig_filtered_early_1.00 <- subset(res_sig_early_1.00,abs(log2FoldChange) > 0)
##res_sig_filtered_early_1.00 <- subset(res_sig_filtered_early_1.00, baseMean > 10)
res_sig_filtered_early_1.00
res_sig_filtered_early_1.00 <- res_sig_filtered_early_1.00[(rownames(res_sig_filtered_early_1.00) %in% differntial_OTU_ID$OTUID),]
res_sig_filtered_early_1.00$Date = "5-21-2018" 
  write.table(res_sig_filtered_late, "/home/hl46161/new_living_mulch/phyloseq_data/differential_OTU_filtered_late", sep="\t")



differntial_OTU_ID_2 <- rbind.data.frame(res_sig_filtered_early_1.00,res_sig_filtered_mid_1.00)
differntial_OTU_ID_2 <-  rbind.data.frame(differntial_OTU_ID_2,res_sig_filtered_late_1.00)
differntial_OTU_ID_2 <- distinct(differntial_OTU_ID)
write.csv(differntial_OTU_ID_2,file = "common_differential_OTU_in_three_time_point_silva.csv",row.names = TRUE,col.names = TRUE,sep = ",")
  ?write.csv


##
increasing_OTU	=	read.csv("/home/hl46161/PycharmProjects/living_mulch/differentially_OTU_ID_increasing.csv",sep=",")
decreasing_OTU	=	read.csv("/home/hl46161/PycharmProjects/living_mulch/differentially_OTU_ID_decreasing.csv",sep=",")
increase_decreasing_OTU	=	read.csv("/home/hl46161/PycharmProjects/living_mulch/differentially_OTU_ID_increas_then_decrease.csv",sep=",")
decreasing_increase_OTU	=	read.csv("/home/hl46161/PycharmProjects/living_mulch/differentially_OTU_ID_decreas_then_increase_OTU.csv",sep=",")

ggplot(data = increasing_OTU,aes(x = Date,y= log2FoldChange, group = OTUID)) + geom_line(aes(color = OTUID)) + geom_point()
ggplot(data = decreasing_OTU,aes(x = Date,y= log2FoldChange, group = OTUID)) + geom_line(aes(color = OTUID)) + geom_point()
ggplot(data = increase_decreasing_OTU,aes(x = Date,y= log2FoldChange, group = OTUID)) + geom_line(aes(color = OTUID)) + geom_point()
ggplot(data = decreasing_increase_OTU,aes(x = Date,y= log2FoldChange, group = OTUID)) + geom_line(aes(color = OTUID)) + geom_point()

str(increasing_OTU)

?rbind.data.frame

##network analysis 

##filtering low variance OTU less than 0.001 change
psr = transform_sample_counts(ps, function(x) x / sum(x) )
psfr = filter_taxa(psr, function(x) mean(x) > 1e-3, TRUE)

#sample_ids <- data.frame(SampleID = row.names(sample_data(ps)))
#tax_table(ps)
#tax_table(psfr)
#row.names(sample_ids) <- sample_ids$SampleID
#exp_sample_ids <- sample_data(sample_ids)
#psfr <- merge_phyloseq(psfr, exp_sample_ids)

#plot_net(psfr, maxdist = 0.4, point_label = "treatment",color = "treatment",shape = "Date_taken")


ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)

##write_graph(OTU_network, file = "OTU_network.txt", "edgelist",)
?write_graph
class(OTU_network)
##network_graph <- plot_net(psfr,type = "taxa",maxdist = 0.3,color = "genus")
##png(filename = "living_mulch_net_analysis.png", width = 2500,height = 2500)
##network_graph
##dev.off()

#subset living_mulch data 
ps.living_mulch <- subset_samples(ps, treatment %in% c("LivingMulch"))
sample_data(ps)
sample_data(ps.living_mulch) 
makeps.living_mulchr = transform_sample_counts(ps.living_mulch, function(x) x / sum(x) )
#filter about low abudance OTU 
ps.living_mulchfr = filter_taxa(ps.living_mulchr, function(x) mean(x) > 1e-3, TRUE)
#ps.living_mulch.rarefied = rarefy_even_depth(ps.living_mulchfr, rngseed=1, sample.size=0.9*min(sample_sums(ps.living_mulch)), replace=F)

#subset No_Cover data
ps.NoCover  <- subset_samples(ps, treatment %in% c("NoCover"))
sample_data(ps)
ps.NoCover
ps.NoCoverr = transform_sample_counts(ps.NoCover, function(x) x / sum(x) )
ps.NoCoverfr = filter_taxa(ps.NoCoverr, function(x) mean(x) > 1e-3, TRUE)


OTU_network <- make_network(psfr,type = "taxa",max.dist = 0.3,color = "genus")
OTU_network
OTU_network_lm <- make_network(ps.living_mulchfr,type = "taxa",max.dist = 0.3,color = "genus")
OTU_network_lm
OTU_network_nc <- make_network(ps.NoCoverfr,type = "taxa",max.dist = 0.3,color = "genus")
OTU_network_nc
createNetworkFromIgraph(OTU_network,"myIgraph.txt")
createNetworkFromIgraph(OTU_network_lm,"Living_mulch_all_time.txt")
createNetworkFromIgraph(OTU_network_nc,"NoCover_all_time.txt")


##use fastsparcc to calculate the covariance and coexpression 
##first need to remove low abundance OTU 
low_OTU = read.csv(file = "/home/hl46161/Documents/low_abundance_OTU.csv",check.names = FALSE)
names(otu_table)
low_OTU_name <- low_OTU$`OTU ID`
otu_table	=	read.csv("living_mulch_filtered_no_Inrow_no_blank.csv",sep="\t",skip =1,row.names =1,check.names = FALSE)
otu_table <- otu_table[!(rownames(otu_table) %in% low_OTU_name),]
?write.csv()
write.table(otu_table, file ="/home/hl46161/new_living_mulch/SparCC/living_mulch_filtered_no_Inrow_no_blank_sparcc.tsv",sep ="\t",)




  ##ggsave("living_mulch_network_analysis",width = 5,height = 5)
?ggsave()
sample_data(ps.rarefied)
sample_names(ps.rarefied)
##Taxonomic analysis plot 
plot_bar(ps.rarefied_abund, x="treatment",fill = "phylum",facet_grid = "Date_taken")+geom_bar(stat="identity")


#linear regression
Shannon_file = read.table("/home/hl46161/G2F_data/core-metrics-results-original-table-2300/alpha-diversity-shannon.tsv",sep ="\t",header = 1)
colnames(Shannon_file)[1] <- "SampleID"
meta_data_file = read.table(file = "/home/hl46161/G2F_data/revised_2018_G2F_metadata.tsv", header=TRUE, sep="\t")
merge_file = merge(Shannon_file,meata_data_file,by = "SampleID")
merge_file <- na.omit(merge_file) 
lm_relation_shannon <- lm(merge_file$shannon ~merge_file$Pedigree,data = merge_file)
summary(lm_relation_shannon)

ggplot(lm_relation_shannon, aes(x = merge_file$treatment, y = merge_file$shannon,color = merge_file$Date_taken)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")
  anova(lm_relation_shannon)
    
##

convert_meta_file <- read.csv("/home/hl46161/new_living_mulch/online_network/living_mulch_data_metafile.csv")
write.table(convert_meta_file ,file = "/home/hl46161/new_living_mulch/online_network/living_mulch_data_metafile.tsv",sep = "\t",row.names = FALSE)
?phyloseq
  phyloseq


