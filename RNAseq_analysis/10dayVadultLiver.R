#cold room experiment liver vs. my experiment liver
#use NYBZ in my data because i think mal only had NYBZ (but at room temp and cold room)
#doing this all in males to make things comparable across groups

#Q1: in terms of regulatory types, is my NYBZ similar to either the cold room or the warm room? I would hope more similar to the cold room but probs not

#Q2: is there variation in fidelity among regulatory types across development? ie. are cis changes more robust across developmental stages? (compare to warm and cold separately)
    #compare 10day NYBZ reg breakdown to NYBZ cold breakdown and NYBZ warm breakdown (this analysis would also address Q1)

#Q3: are genes found to be plastic in the cold room similar to those found to be maternal in my experiment? Goal is to confirm the role of temperature
    #Instead of defining plasticity, could be a correlation between 1) LFC BZ-pure to BZ-CF TO LFC BZ-warm v BZ-cold, 2) LFC NY-pure to NY-CF TO LFC NY-cold v NY-warm


#load libraries 
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(data.table)

#load data
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#10day data
all_countData <- read.csv('all_samples_counts.csv', header = TRUE, sep = ",")
all_metadata <- read.csv('all_samples_metadata.csv', header = TRUE, sep = ",")
all_metadata$id <- gsub("-", ".", all_metadata$id)

#liver only and only males
liver_countDataM_10day <- all_countData[ , colnames(all_countData) %like% "liver" & colnames(all_countData) %like% "M" | colnames(all_countData) %like% "gene" ]
liver_metadataM_10day <- all_metadata %>% filter(tissue == "liver" & sex =="M")

#remove F1s for parent and plasticity only analyses
liver_countDataM_10day_pure <- liver_countDataM_10day[ , colnames(liver_countDataM_10day) %like% "NN" | colnames(liver_countDataM_10day) %like% "BB" | colnames(liver_countDataM_10day) %like% "gene"]
liver_metadataM_10day_pure <- subset(liver_metadataM_10day, pop != "BZNY" & pop != "NYBZ")

#retain only F1s for regulatory analyses
liver_countDataM_10day_F1 <- liver_countDataM_10day[ , colnames(liver_countDataM_10day) %like% "NB" | colnames(liver_countDataM_10day) %like% "BN"]
liver_metadataM_10day_F1 <- subset(liver_metadataM_10day, pop == "BZNY" | pop == "NYBZ")


#read in Mal's data, filter for males
adult_countdata <- read.table('all_parents_counts.txt', header = TRUE, sep = "")
adult_metadata <- read.table('sampleMatrix_liver.txt', header = TRUE, sep = "")
#filter metadata for males - sep for F1 and pureline
adult_metadataM <- adult_metadata %>% filter(sex == "M")
adult_metadataM_pure <- adult_metadata %>% filter(pop != "NYBZ" & sex =="M")
adult_metadataM_F1 <- adult_metadata %>% filter(pop == "NYBZ" & sex =="M")
#filter countdata based on the ID being present in the metadata
col_ids <- sub("\\..*$", "", colnames(adult_countdata))
adult_countdataM_pure <- adult_countdata[, col_ids %in% as.character(adult_metadataM$id) | colnames(adult_countdata) == "gene"] #there is no F1 data in the adult_countdata







######################################################################################################################## #

#####PURELINE PLASTICITY ANALYSES######

######################################################################################################################## #

#10 DAY
##### pure line diffs in 10 day old data
#design builds tests the effects of genotype, maternal environment and their interaction
dds_liver_10day <- DESeqDataSetFromMatrix(countData = liver_countDataM_10day_pure,
                                    colData = liver_metadataM_10day_pure,
                                    design = ~ genotype + maternal + genotype:maternal, tidy = TRUE)
keepAvg_liver_10day <- rowMeans(counts(dds_liver_10day)) >= 10 #filter for a mean of greater than 10 reads
dds_liver_10day <- dds_liver_10day[keepAvg_liver_10day,]

#Defining our comparisons 
#get the model matrix
mod_mat_liver_10day <- model.matrix(design(dds_liver_10day), colData(dds_liver_10day))
#Define coefficient vectors for each condition
#maternal and genotype status separate
NY_pureline <- colMeans(mod_mat_liver_10day[dds_liver_10day$genotype == "NY" & dds_liver_10day$maternal == "NY", ])
BZ_pureline <- colMeans(mod_mat_liver_10day[dds_liver_10day$genotype == "BZ" & dds_liver_10day$maternal == "BZ", ])
NY_cf <- colMeans(mod_mat_liver_10day[dds_liver_10day$genotype == "NY" & dds_liver_10day$maternal == "BZ", ])
BZ_cf <- colMeans(mod_mat_liver_10day[dds_liver_10day$genotype == "BZ" & dds_liver_10day$maternal == "NY", ])



#ADULT
##### pure line diffs in adult data
#design builds tests the effects of genotype, environment and their interaction
dds_liver_adult <- DESeqDataSetFromMatrix(countData = adult_countdataM_pure,
                                          colData = adult_metadataM_pure,
                                          design = ~ pop + Environment + pop:Environment, tidy = TRUE)
keepAvg_liver_adult <- rowMeans(counts(dds_liver_adult)) >= 10 #filter for a mean of greater than 10 reads
dds_liver_adult <- dds_liver_adult[keepAvg_liver_adult,]

#Defining our comparisons 
#get the model matrix
mod_mat_liver_adult <- model.matrix(design(dds_liver_adult), colData(dds_liver_adult))
#Define coefficient vectors for each condition
#maternal and genotype status separate
NY_warm <- colMeans(mod_mat_liver_adult[dds_liver_adult$pop == "NY" & dds_liver_adult$Environment == "Warm", ])
BZ_warm <- colMeans(mod_mat_liver_adult[dds_liver_adult$pop == "BZ" & dds_liver_adult$Environment == "Warm", ])
NY_cold <- colMeans(mod_mat_liver_adult[dds_liver_adult$pop == "NY" & dds_liver_adult$Environment == "Cold", ])
BZ_cold <- colMeans(mod_mat_liver_adult[dds_liver_adult$pop == "BZ" & dds_liver_adult$Environment == "Cold", ])


#run the deseqs
deseq_liver_10day <- DESeq(dds_liver_10day, test = "Wald")
deseq_liver_adult <- DESeq(dds_liver_adult, test = "Wald")


#get CF v non-CF 10 day differences, and warm v cold adult differences

#10 DAY
#Pureline vs cross foster (for NY):
res_NY_pureVcf_liver_10day <- results(deseq_liver_10day, contrast = NY_pureline - NY_cf, pAdjustMethod = "BH") # neg LFC = higher expression NY-CF (warm)
NY_pureVcf_liver_10day <- as.data.frame(res_NY_pureVcf_liver_10day)
NY_pureVcf_liver_10day <- NY_pureVcf_liver_10day %>% dplyr::rename(NYpureVcf.lfc = log2FoldChange, NYpureVcf.padj = padj)

#Pureline vs cross foster (for BZ):
res_BZ_cfVpure_liver_10day <- results(deseq_liver_10day, contrast = BZ_cf - BZ_pureline, pAdjustMethod = "BH") # neg LFC = higher expression BZ-PURELINE!! (warm)
BZ_cfVpure_liver_10day <- as.data.frame(res_BZ_cfVpure_liver_10day)
BZ_cfVpure_liver_10day <- BZ_cfVpure_liver_10day %>% dplyr::rename(BZpcfVpure.lfc = log2FoldChange, BZpcfVpure.padj = padj)

#ADULT
#Pureline vs cross foster (for NY):
res_NY_coldVwarm_liver_adult <- results(deseq_liver_adult, contrast = NY_cold - NY_warm, pAdjustMethod = "BH") # neg LFC = higher expression NY-CF (warm)
NY_coldVwarm_liver_adult <- as.data.frame(res_NY_coldVwarm_liver_adult)
NY_coldVwarm_liver_adult <- NY_coldVwarm_liver_adult %>% dplyr::rename(NYcoldVwarm.lfc = log2FoldChange, NYcoldVwarm.padj = padj)

#Pureline vs cross foster (for BZ):
res_BZ_coldVwarm_liver_adult <- results(deseq_liver_adult, contrast = BZ_cold - BZ_warm, pAdjustMethod = "BH") # neg LFC = higher expression BZ-PURELINE!! (warm)
BZ_coldVwarm_liver_adult <- as.data.frame(res_BZ_coldVwarm_liver_adult)
BZ_coldVwarm_liver_adult <- BZ_coldVwarm_liver_adult %>% dplyr::rename(BZpcoldVwarm.lfc = log2FoldChange, BZpcoldVwarm.padj = padj)



#merge BZ and NY datasets to do a scatterplot of lfcs
liver_mergeNY <- na.omit(merge(NY_pureVcf_liver_10day,NY_coldVwarm_liver_adult, by="row.names"))
liver_mergeBZ <- na.omit(merge(BZ_cfVpure_liver_10day,BZ_coldVwarm_liver_adult, by="row.names"))

#filter merged datasets for only genes that are significantly differentially expressed
liver_mergeNY_sig <- liver_mergeNY %>% filter(NYpureVcf.padj < 0.05 & NYcoldVwarm.padj < 0.05) #whoa only 49 genes (about 2.7% of sig diff NY 10 day genes - but about a 1/4th of cold v warm)
liver_mergeBZ_sig <- liver_mergeBZ %>% filter(BZpcfVpure.padj < 0.05 & BZpcoldVwarm.padj < 0.05) #whoa only 58 genes (about 10% of the sig dif 10 day diffs - but 1/6th of all BZ diffs )
#the differences in expression between cold and warm among BZ and NY are very similar to Mal's numbers

#find separate sig genes
BZ_cfVpure_liver_10day_sig <- BZ_cfVpure_liver_10day %>% filter(BZpcfVpure.padj < 0.05) #595
NY_pureVcf_liver_10day_sig <- NY_pureVcf_liver_10day %>% filter(NYpureVcf.padj < 0.05) #1832

BZ_coldVwarm_liver_adult_sig <- BZ_coldVwarm_liver_adult %>% filter(BZpcoldVwarm.padj < 0.05) #830
NY_coldVwarm_liver_adult_sig <- NY_coldVwarm_liver_adult %>% filter(NYcoldVwarm.padj < 0.05) #204

#ENRICHMENT of plasticity overlap

####enrichmentin sig BZ DEGs
permSigDEGsOverlap.BZliver <- rep(NA, 10000) 
genesToDraw10 <- 595 # number of genes that are sig DEGs BZ cf to non cf 10 day liver comparison
genesToDrawAdult <- 830 # number of genes that are sig DEGs BZ warm to cold adult liver comparison

for(i in 1:length(permSigDEGsOverlap.BZliver)){
  drawA <- sample(x=as.list(row.names(BZ_cfVpure_liver_10day)), size=genesToDraw10, replace=FALSE)
  drawB <- sample(x=as.list(row.names(BZ_coldVwarm_liver_adult)), size=genesToDrawAdult, replace=FALSE)
  permSigDEGsOverlap.BZliver[i] <- length(intersect(drawA, drawB))
}

permSigDEGsOverlap.BZliver <- as.data.frame(permSigDEGsOverlap.BZliver) #we know the overlap is 58
prop(~ permSigDEGsOverlap.BZliver >= 58, data = permSigDEGsOverlap.BZliver) #0 - oh!
#####histogram of simulated distribution
ggplot(permSigDEGsOverlap.BZliver, aes(x=permSigDEGsOverlap.BZliver)) + geom_histogram(color="black", fill="grey89", binwidth=2) + 
  labs(x="Expected overlap", y = "Count") + geom_vline(aes(xintercept=58), color="#ce8c14", linetype="dashed", size=2) +
  theme_classic() +theme(axis.title.x = element_text(face="bold", size =12),
                         axis.text.x = element_text(size =9, colour = 'black'), 
                         axis.title.y = element_text(face="bold", size =12),
                         axis.text.y = element_text(size =9, colour = 'black'))


####enrichmentin sig NY DEGs
permSigDEGsOverlap.NYliver <- rep(NA, 10000) 
genesToDraw10 <- 1832 # number of genes that are sig DEGs NY cf to non cf 10 day liver comparison
genesToDrawAdult <- 204 # number of genes that are sig DEGs NY warm to cold adult liver comparison

for(i in 1:length(permSigDEGsOverlap.NYliver)){
  drawA <- sample(x=as.list(row.names(NY_pureVcf_liver_10day)), size=genesToDraw10, replace=FALSE)
  drawB <- sample(x=as.list(row.names(NY_coldVwarm_liver_adult)), size=genesToDrawAdult, replace=FALSE)
  permSigDEGsOverlap.NYliver[i] <- length(intersect(drawA, drawB))
}

permSigDEGsOverlap.NYliver <- as.data.frame(permSigDEGsOverlap.NYliver) #we know the overlap is 49
prop(~ permSigDEGsOverlap.NYliver >= 49, data = permSigDEGsOverlap.NYliver) #0 - oh!
#####histogram of simulated distribution
ggplot(permSigDEGsOverlap.NYliver, aes(x=permSigDEGsOverlap.NYliver)) + geom_histogram(color="black", fill="grey89", binwidth=2) + 
  labs(x="Expected overlap", y = "Count") + geom_vline(aes(xintercept=49), color="#006d92", linetype="dashed", size=2) +
  theme_classic() +theme(axis.title.x = element_text(face="bold", size =12),
    axis.text.x = element_text(size =9, colour = 'black'), 
    axis.title.y = element_text(face="bold", size =12),
    axis.text.y = element_text(size =9, colour = 'black'))




#SCATTERPLOTS of significant diffs
#correlation between 1) LFC BZ-pure to BZ-CF TO LFC BZ-warm v BZ-cold, 2) LFC NY-pure to NY-CF TO LFC NY-cold v NY-warm

#NY
ggplot(liver_mergeNY_sig, aes(x = NYpureVcf.lfc, y = NYcoldVwarm.lfc)) + xlab("NY 10 day: pure (NY mom-cold) VS. CF (BZ mom-warm)") +
  geom_point(shape=21, fill = "#006d92", color = "black", size=2) +theme_bw() + ylab("NY adult: Cold VS. Warm") +
  geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = -1, label.y = 2.5) #R = 0.02, p = 0.87

#BZ
ggplot(liver_mergeBZ_sig, aes(x = BZpcfVpure.lfc, y = BZpcoldVwarm.lfc)) + xlab("BZ 10 day: CF (NY mom-cold) VS. pure (BZ mom-warm)") +
  geom_point(shape=21, fill = "#ce8c14", color = "black", size=2) +theme_bw() +  ylab("BZ adult: Cold VS. Warm") +
  geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = -1.2, label.y = 2) #R = 0.4, p=0.002 this is sig



#Get normalized counts for each group separately to correlate

#first need a combo deseq option with all of our count data to accurately normalize
adult_10dayCountDataMpure <- na.omit(merge(liver_countDataM_10day_pure,adult_countdataM_pure, by="gene"))
adult_10dayMetadataMpure <- read.csv("adult_10dayMetadataMpure.csv", header = TRUE, sep = ",")

#design builds tests the effects of genotype, environment and their interaction
dds_liver_adult10day <- DESeqDataSetFromMatrix(countData = adult_10dayCountDataMpure,
                                          colData = adult_10dayMetadataMpure,
                                          design = ~ age + genotype + environment + genotype:environment, tidy = TRUE)
keepAvg_liver_adult10day <- rowMeans(counts(dds_liver_adult10day)) >= 10 #filter for a mean of greater than 10 reads
dds_liver_adult10day <- dds_liver_adult10day[keepAvg_liver_adult10day,]

#run the deseq
deseq_liver_adult10day <- DESeq(dds_liver_adult10day, test = "Wald")

#get normalized counts
norm_counts_adult10day <- counts(deseq_liver_adult10day, normalized = TRUE)
groups_adult10day <- colData(deseq_liver_adult10day)[["group"]]
#calculate basemeans per group
base_means_adult10day <- apply(norm_counts_adult10day, 1, function(gene_counts) {tapply(gene_counts, groups_adult10day, mean)})
# Convert to a dataframe
base_means_adult10day <- as.data.frame(t(base_means_adult10day))

#PCA of pureline BZ and NY samples from both experiments
#first normalize and convert to a matrix
vst_all_samples <- vst(dds_liver_adult10day, blind=FALSE)
vst_all_samples_mat <- assay(vst_all_samples)

#make pca object
pca_all_samples <- prcomp(t(vst_all_samples_mat))

#create data frame with metadata and PC1 and PC2 values for input to ggplot
pca_all_samples_df <- as.data.frame(pca_all_samples$x)

#calculate PC percentages
percentage_all <- round((pca_all_samples$sdev^2) / (sum(pca_all_samples$sdev^2)) * 100, 2)
percentage_all <- paste(colnames(pca_all_samples_df), "(",paste(as.character(percentage_all),"%", ")", sep="") )

#add metadata
pca_all_samples_df <- cbind(adult_10dayMetadataMpure, pca_all_samples$x)

#plot
PCA_plot_all <- ggplot(pca_all_samples_df) + geom_point(aes(x=PC1, y=PC2, color=sample_group, shape = age), size =3, alpha=0.9) +
  theme_classic() + xlab(percentage_all[1]) + ylab(percentage_all[2]) +
  theme(legend.position = 'right', 
        axis.line.x = element_line(colour = 'black', size=0.7, linetype='solid'), 
        axis.title.x = element_text(face="bold", size =10),
        axis.text.x = element_text(size =8, colour = 'black'), 
        axis.line.y = element_line(colour = 'black', size=0.7, linetype='solid'),
        axis.title.y = element_text(face="bold", size =10),
        axis.text.y = element_text(size =8, colour = 'black'))
PCA_plot_all

PCA_plot_all_2v3 <- ggplot(pca_all_samples_df) + geom_point(aes(x=PC2, y=PC3, color=sample_group, shape = age), size =3, alpha=0.9) +
  theme_classic() + xlab(percentage_all[2]) + ylab(percentage_all[3]) +
  theme(legend.position = 'right', 
        axis.line.x = element_line(colour = 'black', size=0.7, linetype='solid'), 
        axis.title.x = element_text(face="bold", size =10),
        axis.text.x = element_text(size =8, colour = 'black'), 
        axis.line.y = element_line(colour = 'black', size=0.7, linetype='solid'),
        axis.title.y = element_text(face="bold", size =10),
        axis.text.y = element_text(size =8, colour = 'black'))
PCA_plot_all_2v3

PCA_plot_all_2v4 <- ggplot(pca_all_samples_df) + geom_point(aes(x=PC2, y=PC4, color=sample_group, shape = age), size =3, alpha=0.9) +
  theme_classic() + xlab(percentage_all[2]) + ylab(percentage_all[4]) +
  theme(legend.position = 'right', 
        axis.line.x = element_line(colour = 'black', size=0.7, linetype='solid'), 
        axis.title.x = element_text(face="bold", size =10),
        axis.text.x = element_text(size =8, colour = 'black'), 
        axis.line.y = element_line(colour = 'black', size=0.7, linetype='solid'),
        axis.title.y = element_text(face="bold", size =10),
        axis.text.y = element_text(size =8, colour = 'black'))
PCA_plot_all_2v4


#base mean correlations - these are same in both because basemean is taken from the all adult and all 10day deseqs
#NY - 10day cold v adult cold
ggplot(base_means_adult10day, aes(x = NY_Cold_10day, y = NY_Cold_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#006d92", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.76, very sig
#NY - 10day cold v adult warm (should be less than above)
ggplot(base_means_adult10day, aes(x = NY_Cold_10day, y = NY_Warm_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#006d92", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.77, very sig - not shocking
#BZ - 10day cold v adult cold
ggplot(base_means_adult10day, aes(x = BZ_Cold_10day, y = BZ_Cold_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#ce8c14", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.77, very sig
#BZ - 10day cold v adult warm (should be less than above)
ggplot(base_means_adult10day, aes(x = BZ_Cold_10day, y = BZ_Warm_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#ce8c14", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.78, very sig - not shocking
#BZ v NY - 10day BZ cold to NY adult cold (just seeing if everything is just really correlated)
ggplot(base_means_adult10day, aes(x = BZ_Cold_10day, y = NY_Cold_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#006d92", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.77, its just that liver expression is liver expression
#BZ v NY - 10day BZ warm to NY adult cold (just seeing if everything is just really correlated, because this is the most polar comp)
ggplot(base_means_adult10day, aes(x = BZ_Warm_10day, y = NY_Cold_adult)) + scale_y_continuous(limits = c(0,100000)) +scale_x_continuous(limits = c(0,100000)) + 
  geom_point(shape=21, fill = "#006d92", color = "black") +theme_bw() + geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "solid") +
  stat_cor(method = "pearson", label.x = 1, label.y = 100000) #R = 0.75, its just that liver expression is liver expression
















######################################################################################################################## #

#####REGULATORY ARCHITECTURE ANALYSES#####

######################################################################################################################## #


#NEED TO
#redefine regulatory categories using Mal's data
#and redefine my NYBZ regulatory categories using just males to be comparable

##### 10 day define Parent, Hybrid, LRT tests #####

#read in data - count table that contains each F1 with BZ and NY reads in subsequent columns
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#read in data and metadata
#these files are all samples together, can later be broken up for specific comparisons
F1_ASE_countData <- read.csv('F1_ASE_samples_counts.csv', header = TRUE, sep = ",")
F1_ASE_metadata <- read.csv('F1_ASE_metadata.csv', header = TRUE, sep = ",")
F1_ASE_metadata$id <- gsub("-", ".", F1_ASE_metadata$id)

#liver only 
F1_ASE_liver_countDataM <- F1_ASE_countData[ , colnames(F1_ASE_countData) %like% "liver" & colnames(F1_ASE_countData) %like% "M" | colnames(F1_ASE_countData) %like% "gene"]
F1_ASE_liver_metadataM <- F1_ASE_metadata[F1_ASE_metadata$id %like% "liver" & F1_ASE_metadata$id %like% "M", ]

#separate NYBZ and BZNY
NYBZ_liver_countDataM <- F1_ASE_liver_countDataM[ , colnames(F1_ASE_liver_countDataM) %like% "NB" |  colnames(F1_ASE_liver_countDataM) %like% "gene"]
NYBZ_liver_metadataM <- F1_ASE_liver_metadataM[F1_ASE_liver_metadataM$individual %like% "NB" , ]


###### make dds objects, and run analysis for each separate F1 group

#NYBZ
ddsNYBZ_liver_10day <- DESeqDataSetFromMatrix(countData=NYBZ_liver_countDataM, 
                                  colData=NYBZ_liver_metadataM, 
                                  design=~0+individual+allele, tidy = TRUE)
##set factors size to 1 because we are comparing ratios
sizeFactors(ddsNYBZ_liver_10day) <- rep(1,12) #we have 12 samples
keepNYBZ_liver_10day <- rowMeans(counts(ddsNYBZ_liver_10day)) >= 10
ddsNYBZ_liver_10day <- ddsNYBZ_liver_10day[keepNYBZ_liver_10day,] 
deseq_NYBZ_liver_10day <- DESeq(ddsNYBZ_liver_10day)
res_NYBZ_liver_10day <- results(deseq_NYBZ_liver_10day, contrast=c("allele","NY","BZ"))
NYBZ_ASE_liver_10day <- as.data.frame(res_NYBZ_liver_10day)
#7739 genes to test
####informative names
NYBZ_ASE_liver_10day$gene <- rownames(NYBZ_ASE_liver_10day)
names(NYBZ_ASE_liver_10day)[2]<-paste("NYBZ.NYvBZ_Log2FC")
names(NYBZ_ASE_liver_10day)[6]<-paste("NYBZ.NYvBZ_padj")

####write total ASE genes to tables
write.csv(NYBZ_ASE_liver_10day, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZ_liver10day_ASEgenes.csv", row.names = TRUE)



##### Parental data - test 2: wald's test in DEseq on just parents
#this is analogous to the NYvBZ pureline test above

#Need to 1) subset the original pureline test to only include ASE genes (sep for NYBZ and BZNY) and 
# 2) rename columns to be informative

#NYBZ version
res_parents_NYBZ_liver10day <- results(deseq_liver_10day, contrast = NY_pureline - BZ_pureline, pAdjustMethod = "BH") # neg LFC = higher expression BZ
parents_NYBZ_liver10day <- as.data.frame(res_parents_NYBZ_liver10day)
parents_NYBZ_liver10day$gene <- rownames(parents_NYBZ_liver10day)
parents_NYBZ_liver10day <- subset(parents_NYBZ_liver10day, gene %in% NYBZ_ASE_liver_10day$gene) #7739
names(parents_NYBZ_liver10day)[2]<-paste("parent.NYvBZ_Log2FC")
names(parents_NYBZ_liver10day)[6]<-paste("parent.NYvBZ_padj")



#### Trans component using LRT
#data - test 3: both F1s and parents analyzed together

#LRT based test to see if the population alone (regardless of F1 or parent) explains expression differences. So if it does not, that would point to trans effects being involved
# so if it is significant in this test then the reduced test does not accurately explain differences - which would point to trans, not cis being involved
#use parental counts that are just the genes also found in F1 hyrbids


#Filter parent count data for relevant genes, merge with F1 data - the "C" step takes out CFs
parent_liver_countDataM_10day <- liver_countDataM_10day_pure[ , !colnames(liver_countDataM_10day_pure) %like% "C" | colnames(liver_countDataM_10day_pure) %like% "gene"]

#NYBZ
parent_NYBZ_countDataM_liver10day <- subset(parent_liver_countDataM_10day, gene %in% row.names(NYBZ_ASE_liver_10day))
F1andParent_NYBZ_countDataM_liver10day <- na.omit(merge(parent_NYBZ_countDataM_liver10day, NYBZ_liver_countDataM, by="gene", all=TRUE))

#read in metadata
F1andParent_NYBZ_metadata_liver10day <- read.csv('F1andParent_NYBZ_metadata.csv', header = TRUE, sep = ",")
F1andParent_NYBZ_metadata_liver10day$id <- gsub("-", ".", F1andParent_NYBZ_metadata_liver10day$id)
F1andParent_NYBZ_metadata_liver10day <- F1andParent_NYBZ_metadata_liver10day[F1andParent_NYBZ_metadata_liver10day$id %like% "M" , ]

#### Run LRT test
design = ~F1_parent + population + population:F1_parent

#NYBZ
ddsF1andParent_NYBZ_liver <- DESeqDataSetFromMatrix(F1andParent_NYBZ_countDataM_liver10day, F1andParent_NYBZ_metadata_liver10day, design, tidy = TRUE)
deseq_F1andParent_NYBZ_liver <- DESeq(ddsF1andParent_NYBZ_liver, test="LRT", reduced= ~ population + F1_parent)
resF1andParent_NYBZ_liver <- results(deseq_F1andParent_NYBZ_liver)
#significant here means that the interaction of population:F1_parent has an effect, which in indicative of trans regulation
F1andParent_NYBZ_liver <- as.data.frame(resF1andParent_NYBZ_liver)

#informative names
F1andParent_NYBZ_liver$gene <- rownames(F1andParent_NYBZ_liver)
names(F1andParent_NYBZ_liver)[2]<-paste("int.NYvBZ_Log2FC")
names(F1andParent_NYBZ_liver)[6]<-paste("int.NYvBZ_padj")



#merge the results and LogFC columns of the above test files on gene name
library(plyr)
merge_results_NYBZ_liver10day <- na.omit(join_all(list(parents_NYBZ_liver10day,NYBZ_ASE_liver_10day,F1andParent_NYBZ_liver), by ="gene", type = 'full'))




##### Adult define Parent, Hybrid, LRT tests #####

#once for warm, once for cold

#read in data - count table that contains each F1 with BZ and NY reads in subsequent columns
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#read in data and metadata
#these files are all samples together, can later be broken up for specific comparisons
F1_ASE_countData_adult <- read.table('F1_liver_ASEreads_warmcold_malefemale.txt', header = TRUE, sep = "")
colnames(F1_ASE_countData_adult) <- gsub("X", "", colnames(F1_ASE_countData_adult))
F1_ASE_countDataM_adult <- F1_ASE_countData_adult[ , colnames(F1_ASE_countData_adult) %like% "M" | colnames(F1_ASE_countData_adult) %like% "gene"]
F1_ASE_metadataM_adult <- read.csv('F1_ASE_metadataM_adult.csv', header = TRUE, sep = ",")
F1_ASE_metadataM_adult$individual <- as.character(F1_ASE_metadataM_adult$individual)

#separate warm and cold
NYBZcold_liver_countDataM_adult <- F1_ASE_countDataM_adult[ , colnames(F1_ASE_countDataM_adult) %like% "cold" |  colnames(F1_ASE_countDataM_adult) %like% "gene"]
NYBZcold_liver_metadataM_adult <- F1_ASE_metadataM_adult[F1_ASE_metadataM_adult$environment %like% "cold" , ]
NYBZwarm_liver_countDataM_adult <- F1_ASE_countDataM_adult[ , colnames(F1_ASE_countDataM_adult) %like% "warm" |  colnames(F1_ASE_countDataM_adult) %like% "gene"]
NYBZwarm_liver_metadataM_adult <- F1_ASE_metadataM_adult[F1_ASE_metadataM_adult$environment %like% "warm" , ]


###### make dds objects, and run analysis for each separate F1 group

#NYBZ warm
ddsNYBZwarm <- DESeqDataSetFromMatrix(countData=NYBZwarm_liver_countDataM_adult, 
                                  colData=NYBZwarm_liver_metadataM_adult, 
                                  design=~0+individual+allele, tidy = TRUE)
##set factors size to 1 because we are comparing ratios
sizeFactors(ddsNYBZwarm) <- rep(1,12) #we have 12 samples
keepNYBZwarm <- rowMeans(counts(ddsNYBZwarm)) >= 10
ddsNYBZwarm <- ddsNYBZwarm[keepNYBZwarm,] 
deseq_NYBZwarm <- DESeq(ddsNYBZwarm)
res_NYBZwarm <- results(deseq_NYBZwarm, contrast=c("allele","NY","BZ"))
NYBZwarm_ASE_liver <- as.data.frame(res_NYBZwarm)
#5462 genes to test

#NYBZ cold
ddsNYBZcold <- DESeqDataSetFromMatrix(countData=NYBZcold_liver_countDataM_adult, 
                                  colData=NYBZcold_liver_metadataM_adult, 
                                  design=~0+individual+allele, tidy = TRUE)
##set factors size to 1 because we are comparing ratios
sizeFactors(ddsNYBZcold) <- rep(1,12) #we have 12 samples
keepNYBZcold <- rowMeans(counts(ddsNYBZcold)) >= 10
ddsNYBZcold <- ddsNYBZcold[keepNYBZcold,]
deseq_NYBZcold <- DESeq(ddsNYBZcold)
res_NYBZcold <- results(deseq_NYBZcold, contrast=c("allele","NY","BZ"))
NYBZcold_ASE_liver <- as.data.frame(res_NYBZcold)
#5437 genes to test

####informative names
NYBZwarm_ASE_liver$gene <- rownames(NYBZwarm_ASE_liver)
names(NYBZwarm_ASE_liver)[2]<-paste("NYBZwarm.NYvBZ_Log2FC")
names(NYBZwarm_ASE_liver)[6]<-paste("NYBZwarm.NYvBZ_padj")
NYBZcold_ASE_liver$gene <- rownames(NYBZcold_ASE_liver)
names(NYBZcold_ASE_liver)[2]<-paste("NYBZcold.NYvBZ_Log2FC")
names(NYBZcold_ASE_liver)[6]<-paste("NYBZcold.NYvBZ_padj")

####write total ASE genes to tables
write.csv(NYBZwarm_ASE_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZwarm_liver_ASEgenes.csv", row.names = TRUE)
write.csv(NYBZcold_ASE_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZcold_liver_ASEgenes.csv", row.names = TRUE)


##### Parental data - test 2: wald's test in DEseq on just parents

#Need to 1)test NY v BZ in the cold and warm from DESEQ above
#NYBZ warm version
res_parents_NYBZwarm_liver <- results(deseq_liver_adult, contrast = NY_warm - BZ_warm, pAdjustMethod = "BH") # neg LFC = higher expression BZ
parents_NYBZwarm_liver <- as.data.frame(res_parents_NYBZwarm_liver)
parents_NYBZwarm_liver$gene <- rownames(parents_NYBZwarm_liver)
parents_NYBZwarm_liver <- subset(parents_NYBZwarm_liver, gene %in% NYBZwarm_ASE_liver$gene) #5462
names(parents_NYBZwarm_liver)[2]<-paste("parent.NYvBZwarm_Log2FC")
names(parents_NYBZwarm_liver)[6]<-paste("parent.NYvBZwarm_padj")

#NYBZ cold version
res_parents_NYBZcold_liver <- results(deseq_liver_adult, contrast = NY_cold - BZ_cold, pAdjustMethod = "BH") # neg LFC = higher expression BZ
parents_NYBZcold_liver <- as.data.frame(res_parents_NYBZcold_liver)
parents_NYBZcold_liver$gene <- rownames(parents_NYBZcold_liver)
parents_NYBZcold_liver <- subset(parents_NYBZcold_liver, gene %in% NYBZcold_ASE_liver$gene) #5437
names(parents_NYBZcold_liver)[2]<-paste("parent.NYvBZcold_Log2FC")
names(parents_NYBZcold_liver)[6]<-paste("parent.NYvBZcold_padj")
#numbers of diff exp here similar to Mal's paper


#### Trans component using LRT
#data - test 3: both F1s and parents analyzed together

#Filter parent count data for relevant genes, merge with F1 data - the "C" step takes out CFs
adult_countdataM_warm <- adult_countdata[, col_ids %in% c("MBSD_002","MBSD_003","MBSD_009","MBSD_012","MBSD_013","MBSD_022","MBSD_055","MBSD_061" ,"MBSD_065" ,"MBSD_072","MBSD_075","MBSD_077") | colnames(adult_countdata) == "gene"] 
adult_countdataM_cold <- adult_countdata[, col_ids %in% c("MBSD_005", "MBSD_007", "MBSD_011", "MBSD_015", "MBSD_020" ,"MBSD_024","MBSD_057" ,"MBSD_059" ,"MBSD_063", "MBSD_067" ,"MBSD_069" ,"MBSD_074") | colnames(adult_countdata) == "gene"] 

#NYBZ warm
parent_NYBZwarm_countDataM_liverAdult <- subset(adult_countdataM_warm, gene %in% NYBZwarm_ASE_liver$gene)
F1andParent_NYBZwarm_countDataM_liverAdult <- na.omit(merge(parent_NYBZwarm_countDataM_liverAdult, NYBZwarm_liver_countDataM_adult, by="gene", all=TRUE))

#NYBZ cold
parent_NYBZcold_countDataM_liverAdult <- subset(adult_countdataM_cold, gene %in% NYBZcold_ASE_liver$gene)
F1andParent_NYBZcold_countDataM_liverAdult <- na.omit(merge(parent_NYBZcold_countDataM_liverAdult, NYBZcold_liver_countDataM_adult, by="gene", all=TRUE))

#read in metadata
F1andParent_NYBZwarm_metadata_liverAdult <- read.csv('F1andParent_NYBZwarm_metadata.csv', header = TRUE, sep = ",")
F1andParent_NYBZcold_metadata_liverAdult <- read.csv('F1andParent_NYBZcold_metadata.csv', header = TRUE, sep = ",")
F1andParent_NYBZwarm_metadata_liverAdult$id <- gsub("-", ".", F1andParent_NYBZwarm_metadata_liverAdult$id)
F1andParent_NYBZcold_metadata_liverAdult$id <- gsub("-", ".", F1andParent_NYBZcold_metadata_liverAdult$id)

#### Run LRT test
design = ~F1_parent + population + population:F1_parent

#NYBZ warm 
ddsF1andParent_NYBZwarm_liver <- DESeqDataSetFromMatrix(F1andParent_NYBZwarm_countDataM_liverAdult, F1andParent_NYBZwarm_metadata_liverAdult, design, tidy = TRUE)
deseq_F1andParent_NYBZwarm_liver <- DESeq(ddsF1andParent_NYBZwarm_liver, test="LRT", reduced= ~ population + F1_parent)
resF1andParent_NYBZwarm_liver <- results(deseq_F1andParent_NYBZwarm_liver)
#significant here means that the interaction of population:F1_parent has an effect, which in indicative of trans regulation
F1andParent_NYBZwarm_liver <- as.data.frame(resF1andParent_NYBZwarm_liver)
#informative names
F1andParent_NYBZwarm_liver$gene <- rownames(F1andParent_NYBZwarm_liver)
names(F1andParent_NYBZwarm_liver)[2]<-paste("int.warmNYvBZ_Log2FC")
names(F1andParent_NYBZwarm_liver)[6]<-paste("int.warmNYvBZ_padj")

#NYBZ cold
ddsF1andParent_NYBZcold_liver <- DESeqDataSetFromMatrix(F1andParent_NYBZcold_countDataM_liverAdult, F1andParent_NYBZcold_metadata_liverAdult, design, tidy = TRUE)
deseq_F1andParent_NYBZcold_liver <- DESeq(ddsF1andParent_NYBZcold_liver, test="LRT", reduced= ~ population + F1_parent)
resF1andParent_NYBZcold_liver <- results(deseq_F1andParent_NYBZcold_liver)
#significant here means that the interaction of population:F1_parent has an effect, which in indicative of trans regulation
F1andParent_NYBZcold_liver <- as.data.frame(resF1andParent_NYBZcold_liver)
#informative names
F1andParent_NYBZcold_liver$gene <- rownames(F1andParent_NYBZcold_liver)
names(F1andParent_NYBZcold_liver)[2]<-paste("int.coldNYvBZ_Log2FC")
names(F1andParent_NYBZcold_liver)[6]<-paste("int.coldNYvBZ_padj")


#merge the results and LogFC columns of the above test files on gene name
library(plyr)
merge_results_NYBZwarm_liverAdult <- na.omit(join_all(list(parents_NYBZwarm_liver,NYBZwarm_ASE_liver,F1andParent_NYBZwarm_liver), by ="gene", type = 'full'))
merge_results_NYBZcold_liverAdult <- na.omit(join_all(list(parents_NYBZcold_liver,NYBZcold_ASE_liver,F1andParent_NYBZcold_liver), by ="gene", type = 'full'))

#





##### 10 day define cis and trans #####

#only doing NYBZ because that is what Mal has

#main dataset is "merge_results_NYBZ_liver10day"

#cis only: sig in P and H, not sig in T
#trans only: sig in P, not sig in H, sig in T
#cis + trans: sig in P, H and T. sign of log2FC in P and H is the SAME direction, same direction
#cis x trans: sig in P, H and T. sign of log2FC in P and H is the OPPOSITE direction, opposing 
#compensatory: sig in H, not sig in P, sig in T. no expression divergence due to compensation
#conserved: not sig in H, P, or T. no expression divergence
#ambiguous: the left over potential patterns.

##### CIS ONLY: sig in H and P, not sig in T
cisOnly_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj > 0.05))

##### TRANS ONLY: sig in P, not sig in H, sig in T
transOnly_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj < 0.05))

##### CIS + TRANS: sig in P, H and T. sign of log2FC in P and H is the SAME direction, same direction 
cisANDtrans_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj < 0.05) &
                                   (((merge_results_NYBZ_liver10day$parent.NYvBZ_Log2FC < 0) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_Log2FC < 0)) | ((merge_results_NYBZ_liver10day$parent.NYvBZ_Log2FC > 0) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_Log2FC > 0))))

##### CIS x TRANS: sig in P, H and T. sign of log2FC in P and H is the OPPOSITE direction, opposing 
cisXtrans_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj < 0.05) &
                                 (((merge_results_NYBZ_liver10day$parent.NYvBZ_Log2FC < 0) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_Log2FC > 0)) | ((merge_results_NYBZ_liver10day$parent.NYvBZ_Log2FC > 0) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_Log2FC < 0))))

####### COMPENSATORY: sig in H, not sig in P, sig in T. no expression divergence due to compensation 
compensatory_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj < 0.05))

##### CONSERVED: not sig in H, P, or T. no expression divergence 
conserved_NYBZ_liver10day <- subset(merge_results_NYBZ_liver10day, (merge_results_NYBZ_liver10day$parent.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver10day$NYBZ.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver10day$int.NYvBZ_padj > 0.05))


#### make combo dataset 
#has genes with their mode of regulation specified 
cisOnly_NYBZ_liver10day$regulation_mode <- "cis only"
transOnly_NYBZ_liver10day$regulation_mode <- "trans only"
cisANDtrans_NYBZ_liver10day$regulation_mode <- "cis + trans"
cisXtrans_NYBZ_liver10day$regulation_mode <- "cis x trans"
compensatory_NYBZ_liver10day$regulation_mode <- "compensatory"
conserved_NYBZ_liver10day$regulation_mode <- "conserved"

#bind the above datasets to have one for each comparison (NYvBZ.liver, NYvBZ.BAT. CANvBZ.liver, CANvBZ.BAT)
NYBZ.liver10day_RegulationModes <- rbind(cisOnly_NYBZ_liver10day, transOnly_NYBZ_liver10day, cisANDtrans_NYBZ_liver10day, cisXtrans_NYBZ_liver10day, compensatory_NYBZ_liver10day, conserved_NYBZ_liver10day) 

#write these datasets to my computer for later access 
write.csv(NYBZ.liver10day_RegulationModes, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZ.liver10day_RegulationModes.csv", row.names = TRUE)






##### Adult define cis and trans #####

#once for warm and once for cold

#main datasets are "merge_results_NYBZwarm_liverAdult" and "merge_results_NYBZcold_liverAdult"

##### CIS ONLY: sig in H and P, not sig in T
cisOnly_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj > 0.05))
cisOnly_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj < 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj < 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj > 0.05))

##### TRANS ONLY: sig in P, not sig in H, sig in T
transOnly_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj > 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj < 0.05))
transOnly_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj < 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj > 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj < 0.05))

##### CIS + TRANS: sig in P, H and T. sign of log2FC in P and H is the SAME direction, same direction 
cisANDtrans_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj < 0.05) &
                                        (((merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_Log2FC < 0) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_Log2FC < 0)) | ((merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_Log2FC > 0) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_Log2FC > 0))))
cisANDtrans_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj < 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj < 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj < 0.05) &
                                            (((merge_results_NYBZcold_liverAdult$parent.NYvBZcold_Log2FC < 0) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_Log2FC < 0)) | ((merge_results_NYBZcold_liverAdult$parent.NYvBZcold_Log2FC > 0) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_Log2FC > 0))))

##### CIS x TRANS: sig in P, H and T. sign of log2FC in P and H is the OPPOSITE direction, opposing 
cisXtrans_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj < 0.05) &
                                      (((merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_Log2FC < 0) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_Log2FC > 0)) | ((merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_Log2FC > 0) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_Log2FC < 0))))
cisXtrans_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj < 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj < 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj < 0.05) &
                                          (((merge_results_NYBZcold_liverAdult$parent.NYvBZcold_Log2FC < 0) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_Log2FC > 0)) | ((merge_results_NYBZcold_liverAdult$parent.NYvBZcold_Log2FC > 0) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_Log2FC < 0))))

####### COMPENSATORY: sig in H, not sig in P, sig in T. no expression divergence due to compensation 
compensatory_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj > 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj < 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj < 0.05))
compensatory_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj > 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj < 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj < 0.05))

##### CONSERVED: not sig in H, P, or T. no expression divergence 
conserved_NYBZwarm_liverAdult <- subset(merge_results_NYBZwarm_liverAdult, (merge_results_NYBZwarm_liverAdult$parent.NYvBZwarm_padj > 0.05) & (merge_results_NYBZwarm_liverAdult$NYBZwarm.NYvBZ_padj > 0.05) & (merge_results_NYBZwarm_liverAdult$int.warmNYvBZ_padj > 0.05))
conserved_NYBZcold_liverAdult <- subset(merge_results_NYBZcold_liverAdult, (merge_results_NYBZcold_liverAdult$parent.NYvBZcold_padj > 0.05) & (merge_results_NYBZcold_liverAdult$NYBZcold.NYvBZ_padj > 0.05) & (merge_results_NYBZcold_liverAdult$int.coldNYvBZ_padj > 0.05))


#### make combo dataset 
#has genes with their mode of regulation specified 
cisOnly_NYBZwarm_liverAdult$regulation_mode <- "cis only"
transOnly_NYBZwarm_liverAdult$regulation_mode <- "trans only"
cisANDtrans_NYBZwarm_liverAdult$regulation_mode <- "cis + trans"
cisXtrans_NYBZwarm_liverAdult$regulation_mode <- "cis x trans"
compensatory_NYBZwarm_liverAdult$regulation_mode <- "compensatory"
conserved_NYBZwarm_liverAdult$regulation_mode <- "conserved"

cisOnly_NYBZcold_liverAdult$regulation_mode <- "cis only"
transOnly_NYBZcold_liverAdult$regulation_mode <- "trans only"
cisANDtrans_NYBZcold_liverAdult$regulation_mode <- "cis + trans"
cisXtrans_NYBZcold_liverAdult$regulation_mode <- "cis x trans"
compensatory_NYBZcold_liverAdult$regulation_mode <- "compensatory"
conserved_NYBZcold_liverAdult$regulation_mode <- "conserved"

#bind the above datasets to have one for each comparison (NYvBZ.liver, NYvBZ.BAT. CANvBZ.liver, CANvBZ.BAT)
NYBZwarm.liverAdult_RegulationModes <- rbind(cisOnly_NYBZwarm_liverAdult, transOnly_NYBZwarm_liverAdult, cisANDtrans_NYBZwarm_liverAdult, cisXtrans_NYBZwarm_liverAdult, compensatory_NYBZwarm_liverAdult, conserved_NYBZwarm_liverAdult) 
NYBZcold.liverAdult_RegulationModes <- rbind(cisOnly_NYBZcold_liverAdult, transOnly_NYBZcold_liverAdult, cisANDtrans_NYBZcold_liverAdult, cisXtrans_NYBZcold_liverAdult, compensatory_NYBZcold_liverAdult, conserved_NYBZcold_liverAdult) 

#write these datasets to my computer for later access 
write.csv(NYBZwarm.liverAdult_RegulationModes, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZwarm.liverAdult_RegulationModes.csv", row.names = TRUE)
write.csv(NYBZcold.liverAdult_RegulationModes, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZcold.liverAdult_RegulationModes.csv", row.names = TRUE)




######## 10 day to Adult Regulation Alluvials #######

#Q1: in terms of regulatory types, is my NYBZ similar to either the cold room or the warm room? I would hope more similar to the cold room but probs not

#Q2: is there variation in fidelity among regulatory types across development? ie. are cis changes more robust across developmental stages? (compare to warm and cold separately)
#compare 10day NYBZ reg breakdown to NYBZ cold breakdown and NYBZ warm breakdown (this analysis would also address Q1)



#10 day to Adult cold

#all genes in regulatory data - because there is less agreement across datasets here, I am going to filter out NA versus changing them to conserved
Adult10day_merge_cold <- dplyr::full_join(NYBZ.liver10day_RegulationModes,NYBZcold.liverAdult_RegulationModes, by = "gene") %>%
  select("gene", "category_NYBZ10day" = "regulation_mode.x", "category_NYBZcold" = "regulation_mode.y") %>% #just get the columns we care about
  drop_na() #1473 with NA removed

#add a column that contains the regulatory type for both hybrid comparisons
Adult10day_merge_cold$both <- paste(Adult10day_merge_cold$category_NYBZ10day, Adult10day_merge_cold$category_NYBZcold, sep="_")

#summarize data in plottable format
Adult10day_merge_cold_formatted <- Adult10day_merge_cold %>%
  dplyr::group_by(category_NYBZ10day, category_NYBZcold, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_NYBZ10day = factor(category_NYBZ10day, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#Alluvial plot of fidelity of regulatory types across development - NYBZ 10 day vs NYBZ cold room
Adult10day_cold_alluvial <- Adult10day_merge_cold_formatted %>%
  ggplot(aes(y = n, axis1 = category_NYBZ10day, axis2 = category_NYBZcold)) +
  geom_alluvium(aes(fill= category_NYBZ10day), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.95, width = 0.33, color = "transparent") + ylab("Number of genes") +
  theme_classic() + scale_x_continuous(breaks = c(1, 2), labels = c("SARAxMANA\n10 day", "SARAxMANA\n Adult - Cold")) +
  theme(legend.position = 'none',
        axis.text.x = element_text(face="bold", size =9, colour = 'black'), 
        axis.title.y = element_text(face="bold", size =11),
        axis.text.y = element_text(size =8, colour = 'black'))
Adult10day_cold_alluvial



#10day to Adult warm
#all genes in regulatory data - because there is less agreement across datasets here, I am going to filter out NA versus changing them to conserved
Adult10day_merge_warm <- dplyr::full_join(NYBZ.liver10day_RegulationModes,NYBZwarm.liverAdult_RegulationModes, by = "gene") %>%
  select("gene", "category_NYBZ10day" = "regulation_mode.x", "category_NYBZwarm" = "regulation_mode.y") %>% #just get the columns we care about
  drop_na() #1473 with NA removed

#add a column that contains the regulatory type for both hybrid comparisons
Adult10day_merge_warm$both <- paste(Adult10day_merge_warm$category_NYBZ10day, Adult10day_merge_warm$category_NYBZwarm, sep="_")

#summarize data in plottable format
Adult10day_merge_warm_formatted <- Adult10day_merge_warm %>%
  dplyr::group_by(category_NYBZ10day, category_NYBZwarm, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_NYBZ10day = factor(category_NYBZ10day, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#Alluvial plot of fidelity of regulatory types across development - NYBZ 10 day vs NYBZ warm room
Adult10day_warm_alluvial <- Adult10day_merge_warm_formatted %>%
  ggplot(aes(y = n, axis1 = category_NYBZ10day, axis2 = category_NYBZwarm)) +
  geom_alluvium(aes(fill= category_NYBZ10day), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.95, width = 0.33, color = "transparent") + ylab("Number of genes") +
  theme_classic() + scale_x_continuous(breaks = c(1, 2), labels = c("SARAxMANA\n10 day", "SARAxMANA\n Adult")) +
  theme(legend.position = 'none',
        axis.text.x = element_text(face="bold", size =9, colour = 'black'), 
        axis.title.y = element_text(face="bold", size =11),
        axis.text.y = element_text(size =8, colour = 'black'))
Adult10day_warm_alluvial



#Alluvials of just genes that are cis only or trans only in both
#acutally that isnt informative because most just goes to conserved and you can clearly see the big cis block in both







































