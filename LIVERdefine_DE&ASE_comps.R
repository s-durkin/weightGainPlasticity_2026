
############################################################## #
#### Load and prepare data ####
############################################################## #

#load packages
library(DESeq2)
library(ggplot2)
library(data.table)
library(dplyr)


#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#read in data and metadata
#these files are all samples together, can later be broken up for specific comparisons
all_countData <- read.csv('all_samples_counts.csv', header = TRUE, sep = ",")
all_metadata <- read.csv('all_samples_metadata.csv', header = TRUE, sep = ",")
all_metadata$id <- gsub("-", ".", all_metadata$id)

#liver only 
liver_countData <- all_countData[ , colnames(all_countData) %like% "liver" | colnames(all_countData) %like% "gene"]
liver_metadata <- subset(all_metadata, tissue == "liver")

#remove F1s - these first analyses are just concerning CF vs non-CF
liver_countData <- liver_countData[ , colnames(liver_countData) %like% "NN" | colnames(liver_countData) %like% "BB" | colnames(liver_countData) %like% "gene"]
liver_metadata <- subset(liver_metadata, pop != "BZNY" & pop != "NYBZ")



############################################################## #
#### Make DDS objects ####
############################################################## #

#design builds tests the effects of genotype, maternal environment and their interaction
dds_liver <- DESeqDataSetFromMatrix(countData = liver_countData,
                                   colData = liver_metadata,
                                   design = ~ genotype + maternal + genotype:maternal, tidy = TRUE)
keepAvg_liver <- rowMeans(counts(dds_liver)) >= 10 #filter for a mean of greater than 10 reads
dds_liver <- dds_liver[keepAvg_liver,]



############################################################## #
#### Defining our comparisons ####
############################################################## #

#get the model matrix
mod_mat_liver <- model.matrix(design(dds_liver), colData(dds_liver))

#Define coefficient vectors for each condition
#- the estimated effects of different experimental conditions on gene expression, as determined by the model
#one for each group that could be compared, this is just the column means for all the samples that fit this parameter

#maternal and genotype status separate
NY_pureline <- colMeans(mod_mat_liver[dds_liver$genotype == "NY" & dds_liver$maternal == "NY", ])
BZ_pureline <- colMeans(mod_mat_liver[dds_liver$genotype == "BZ" & dds_liver$maternal == "BZ", ])
NY_cf <- colMeans(mod_mat_liver[dds_liver$genotype == "NY" & dds_liver$maternal == "BZ", ])
BZ_cf <- colMeans(mod_mat_liver[dds_liver$genotype == "BZ" & dds_liver$maternal == "NY", ])

#Genotype sets
NY_geno <- colMeans(mod_mat_liver[dds_liver$genotype == "NY", ]) #will include NY reared by NY and BZ
BZ_geno <- colMeans(mod_mat_liver[dds_liver$genotype == "BZ", ]) #will include BZ reared by NY and BZ

#Maternal sets
NY_maternal <- colMeans(mod_mat_liver[dds_liver$maternal == "NY", ]) #will include NY and BZ reared by NY 
BZ_maternal <- colMeans(mod_mat_liver[dds_liver$maternal == "BZ", ]) #will include NY and BZ reared by BZ 



############################################################## #
#### Run DESeq ####
############################################################## #

deseq_liver <- DESeq(dds_liver, test = "Wald")



############################################################## #
#### Define results summaries ####
############################################################## #

#New York vs Brazil (pureline):
res_pureline_NYvsBZ_liver <- results(deseq_liver, contrast = NY_pureline - BZ_pureline, pAdjustMethod = "BH") # neg LFC = higher expression in BZ 
pureline_NYvsBZ_liver <- as.data.frame(res_pureline_NYvsBZ_liver)

#New York vs Brazil (cross foster):
res_cf_NYvsBZ_liver <- results(deseq_liver, contrast = NY_cf - BZ_cf, pAdjustMethod = "BH") # neg LFC = higher expression in BZ 
cf_NYvsBZ_liver <- as.data.frame(res_cf_NYvsBZ_liver)

#Pureline vs cross foster (for NY):
res_NY_pureVcf_liver <- results(deseq_liver, contrast = NY_pureline - NY_cf, pAdjustMethod = "BH") # neg LFC = higher expression NY-CF
NY_pureVcf_liver <- as.data.frame(res_NY_pureVcf_liver)

#Pureline vs cross foster (for BZ):
res_BZ_cfVpure_liver <- results(deseq_liver, contrast = BZ_cf - BZ_pureline, pAdjustMethod = "BH") # neg LFC = higher expression BZ-PURELINE!!
BZ_cfVpure_liver <- as.data.frame(res_BZ_cfVpure_liver)

#BZ vs NY-CF (BZ mom variation):
res_NYCFvsBZ_liver <- results(deseq_liver, contrast = NY_cf - BZ_pureline, pAdjustMethod = "BH") # neg LFC = higher expression BZ-PURELINE!!
NYCFvsBZ_liver <- as.data.frame(res_NYCFvsBZ_liver)

#NY vs BZ-CF (NY mom variation):
res_NYvsBZCF_liver <- results(deseq_liver, contrast = NY_pureline - BZ_cf, pAdjustMethod = "BH") # neg LFC = higher expression BZ-CF!!
NYvsBZCF_liver <- as.data.frame(res_NYvsBZCF_liver)

#Effect of genotype ("G"):
res_G_liver <- results(deseq_liver, contrast = NY_geno - BZ_geno, pAdjustMethod = "BH") # neg LFC = higher expression in BZ 
liver_G <- as.data.frame(res_G_liver)

#Effect of Maternal Environment ("E"):
res_E_liver <- results(deseq_liver, contrast = NY_maternal - BZ_maternal, pAdjustMethod = "BH") # neg LFC = higher expression in BZ reared mice
liver_E <- as.data.frame(res_E_liver)

#DESEQ's definition of genotype-by-environment ("GxE")
#difference between NYvBZ with NY mom and NYvBZ with BZ mom
res_GxE_liver <- results(deseq_liver, contrast = (NY_pureline - BZ_cf) - (NY_cf - BZ_pureline), pAdjustMethod = "BH")
liver_GxE <- as.data.frame(res_GxE_liver)



#### save R data ####
#this is just the DE comps, no ASE yet

#save the R environment so this can be loaded for all downstream analyses

save.image(file = "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/liver_DE_comps.RData")











############################################################## #
#### ASE ANALYSES ####
############################################################## #



############################################################## #
#### Hybrid alone test ####
############################################################## #

### F1 data - test 1: wald's test in DESeq for differential expression between BZ and NY alleles 

#read in data - count table that contains each F1 with BZ and NY reads in subsequent columns
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#read in data and metadata
#these files are all samples together, can later be broken up for specific comparisons
F1_ASE_countData <- read.csv('F1_ASE_samples_counts.csv', header = TRUE, sep = ",")
F1_ASE_metadata <- read.csv('F1_ASE_metadata.csv', header = TRUE, sep = ",")
F1_ASE_metadata$id <- gsub("-", ".", F1_ASE_metadata$id)

#liver only 
F1_ASE_liver_countData <- F1_ASE_countData[ , colnames(F1_ASE_countData) %like% "liver" | colnames(F1_ASE_countData) %like% "gene"]
F1_ASE_liver_metadata <- subset(F1_ASE_metadata, tissue == "liver")

##separate NYBZ and BZNY
NYBZ_liver_countData <- F1_ASE_liver_countData[ , colnames(F1_ASE_liver_countData) %like% "NB" |  colnames(F1_ASE_liver_countData) %like% "gene"]
NYBZ_liver_metadata <- F1_ASE_liver_metadata[F1_ASE_liver_metadata$individual %like% "NB" , ]

BZNY_liver_countData <- F1_ASE_liver_countData[ , colnames(F1_ASE_liver_countData) %like% "BN" |  colnames(F1_ASE_liver_countData) %like% "gene"]
BZNY_liver_metadata <- F1_ASE_liver_metadata[F1_ASE_liver_metadata$individual %like% "BN" , ]

  

###### make dds objects, and run analysis for each separate F1 group
#NYBZ
ddsNYBZ_liver <- DESeqDataSetFromMatrix(countData=NYBZ_liver_countData, 
                                  colData=NYBZ_liver_metadata, 
                                  design=~0+individual+allele, tidy = TRUE)
##set factors size to 1 because we are comparing ratios
sizeFactors(ddsNYBZ_liver) <- rep(1,24) #we have 24 samples
keepNYBZ_liver <- rowMeans(counts(ddsNYBZ_liver)) >= 10
ddsNYBZ_liver <- ddsNYBZ_liver[keepNYBZ_liver,] 
deseq_NYBZ_liver <- DESeq(ddsNYBZ_liver)
res_NYBZ_liver <- results(deseq_NYBZ_liver, contrast=c("allele","NY","BZ"))
NYBZ_ASE_liver <- as.data.frame(res_NYBZ_liver)
#7679 genes to test

#BZNY
ddsBZNY_liver <- DESeqDataSetFromMatrix(countData=BZNY_liver_countData, 
                                  colData=BZNY_liver_metadata, 
                                  design=~0+individual+allele, tidy = TRUE)
##set factors size to 1 because we are comparing ratios
sizeFactors(ddsBZNY_liver) <- rep(1,24) #we have 24 samples
keepBZNY_liver <- rowMeans(counts(ddsBZNY_liver)) >= 10
ddsBZNY_liver <- ddsBZNY_liver[keepBZNY_liver,]
deseq_BZNY_liver <- DESeq(ddsBZNY_liver)
res_BZNY_liver <- results(deseq_BZNY_liver, contrast=c("allele","NY","BZ"))
BZNY_ASE_liver <- as.data.frame(res_BZNY_liver)
#7471 genes to test


####informative names
NYBZ_ASE_liver$gene <- rownames(NYBZ_ASE_liver)
names(NYBZ_ASE_liver)[2]<-paste("NYBZ.NYvBZ_Log2FC")
names(NYBZ_ASE_liver)[6]<-paste("NYBZ.NYvBZ_padj")
BZNY_ASE_liver$gene <- rownames(BZNY_ASE_liver)
names(BZNY_ASE_liver)[2]<-paste("BZNY.NYvBZ_Log2FC")
names(BZNY_ASE_liver)[6]<-paste("BZNY.NYvBZ_padj")

####write total ASE genes to tables
write.csv(NYBZ_ASE_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZ_liver_ASEgenes.csv", row.names = TRUE)
write.csv(BZNY_ASE_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/BZNY_liver_ASEgenes.csv", row.names = TRUE)


############################################################## #
#### Parent alone test ####
############################################################## #

##### Parental data - test 2: wald's test in DEseq on just parents
#this is analogous to the NYvBZ pureline test above

#Need to 1) subset the original pureline test to only include ASE genes (sep for NYBZ and BZNY) and 
# 2) rename columns to be informative

#NYBZ version
parents_NYBZ_liver <- pureline_NYvsBZ_liver
parents_NYBZ_liver$gene <- rownames(parents_NYBZ_liver)
parents_NYBZ_liver <- subset(parents_NYBZ_liver, gene %in% NYBZ_ASE_liver$gene) #7679

names(parents_NYBZ_liver)[2]<-paste("parent.NYvBZ_Log2FC")
names(parents_NYBZ_liver)[6]<-paste("parent.NYvBZ_padj")

#BZNY version
parents_BZNY_liver <- pureline_NYvsBZ_liver
parents_BZNY_liver$gene <- rownames(parents_BZNY_liver)
parents_BZNY_liver <- subset(parents_BZNY_liver, gene %in% BZNY_ASE_liver$gene) #7471

names(parents_BZNY_liver)[2]<-paste("parent.NYvBZ_Log2FC")
names(parents_BZNY_liver)[6]<-paste("parent.NYvBZ_padj")



############################################################## #
#### Trans component using LRT ####
############################################################## #

#data - test 3: both F1s and parents analyzed together

#LRT based test to see if the population alone (regardless of F1 or parent) explains expression differences. So if it does not, that would point to trans effects being involved
# so if it is significant in this test then the reduced test does not accurately explain differences - which would point to trans, not cis being involved
#use parental counts that are just the genes also found in F1 hyrbids


#Filter parent count data for relevant genes, merge with F1 data
parent_liver_countData <- liver_countData[ , !colnames(liver_countData) %like% "C" | colnames(liver_countData) %like% "gene"]

#NYBZ
parent_NYBZ_countData_liver <- subset(parent_liver_countData, gene %in% NYBZ_ASE_liver$gene)
F1andParent_NYBZ_countData_liver <- merge(parent_NYBZ_countData_liver, NYBZ_liver_countData, by="gene", all=TRUE)
F1andParent_NYBZ_countData_liver <- na.omit(F1andParent_NYBZ_countData_liver)

#BZNY
parent_BZNY_countData_liver <- subset(parent_liver_countData, gene %in% BZNY_ASE_liver$gene)
F1andParent_BZNY_countData_liver <- merge(parent_BZNY_countData_liver, BZNY_liver_countData, by="gene", all=TRUE)
F1andParent_BZNY_countData_liver <- na.omit(F1andParent_BZNY_countData_liver)

#read in metadata
F1andParent_NYBZ_metadata_liver <- read.csv('F1andParent_NYBZ_metadata.csv', header = TRUE, sep = ",")
F1andParent_BZNY_metadata_liver <- read.csv('F1andParent_BZNY_metadata.csv', header = TRUE, sep = ",")
F1andParent_NYBZ_metadata_liver$id <- gsub("-", ".", F1andParent_NYBZ_metadata_liver$id)
F1andParent_BZNY_metadata_liver$id <- gsub("-", ".", F1andParent_BZNY_metadata_liver$id)

#### Run LRT test

design = ~F1_parent + population + population:F1_parent

#NYBZ
ddsF1andParent_NYBZ_liver <- DESeqDataSetFromMatrix(F1andParent_NYBZ_countData_liver, F1andParent_NYBZ_metadata_liver, design, tidy = TRUE)
deseq_F1andParent_NYBZ_liver <- DESeq(ddsF1andParent_NYBZ_liver, test="LRT", reduced= ~ population + F1_parent)
resF1andParent_NYBZ_liver <- results(deseq_F1andParent_NYBZ_liver)
#significant here means that the interaction of population:F1_parent has an effect, which in indicative of trans regulation
F1andParent_NYBZ_liver <- as.data.frame(resF1andParent_NYBZ_liver)

#informative names
F1andParent_NYBZ_liver$gene <- rownames(F1andParent_NYBZ_liver)
names(F1andParent_NYBZ_liver)[2]<-paste("int.NYvBZ_Log2FC")
names(F1andParent_NYBZ_liver)[6]<-paste("int.NYvBZ_padj")


#BZNY
ddsF1andParent_BZNY_liver <- DESeqDataSetFromMatrix(F1andParent_BZNY_countData_liver, F1andParent_BZNY_metadata_liver, design, tidy = TRUE)
deseq_F1andParent_BZNY_liver <- DESeq(ddsF1andParent_BZNY_liver, test="LRT", reduced= ~ population + F1_parent)
resF1andParent_BZNY_liver <- results(deseq_F1andParent_BZNY_liver)
#significant here means that the interaction of population:F1_parent has an effect, which in indicative of trans regulation
F1andParent_BZNY_liver <- as.data.frame(resF1andParent_BZNY_liver)

#informative names
F1andParent_BZNY_liver$gene <- rownames(F1andParent_BZNY_liver)
names(F1andParent_BZNY_liver)[2]<-paste("int.NYvBZ_Log2FC")
names(F1andParent_BZNY_liver)[6]<-paste("int.NYvBZ_padj")





############################################################## #
#### Merge and save final data frames ####
############################################################## #


#result dataframe names for each test type
#P: parent alone wald's test: parents_NYBZ_liver, parents_BZNY_liver
#H: hybrid alone wald's test: NYBZ_ASE_liver, BZNY_ASE_liver
#T: trans component LRT test: F1andParent_NYBZ_liver, F1andParent_BZNY_liver

#merge the results and LogFC columns of the above test files on gene name
library(plyr)
merge_results_NYBZ_liver <- join_all(list(parents_NYBZ_liver,NYBZ_ASE_liver,F1andParent_NYBZ_liver), by ="gene", type = 'full')
merge_results_NYBZ_liver <- na.omit(merge_results_NYBZ_liver) #7679 genes

merge_results_BZNY_liver <- join_all(list(parents_BZNY_liver,BZNY_ASE_liver,F1andParent_BZNY_liver), by ="gene", type = 'full')
merge_results_BZNY_liver <- na.omit(merge_results_BZNY_liver) #7471 genes




#### save R data ####
#this includes DE and ASE comparisons

#save the R environment so this can be loaded for all downstream cis v trans analyses

save.image(file = "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/liver_ASE_comps.RData")


