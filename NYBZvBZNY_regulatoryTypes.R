#NYBZ vs BZNY regulation across gene sets

#the idea is that NYBZ is a "cold environment" hybrid, and BZNY the warm one. 

#Intereseting to think about how regulatory categories shift between the two hybrids - as a function of temperature

#Does the stability of cis vs trans changes differ across tissues and across gene sets (genetic vs. maternal vs. all genes)

#Plots similar to Mal's PNAS paper Figure 3C
#Her github contains helpful code here: https://github.com/malballinger/BallingerMack_PNAS_2023/blob/main/code/postprocess_RNAseq/plot_cistrans.R

#Load packages
library(dplyr)
library(ggplot2)
library(tidyverse)

#read in Regulatory data for NYBZ and BZNY and the genetic and maternal genes
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#load regulatory modes spreadsheets - both tissues
NYBZ.tail_RegulationModes <- read.csv("NYBZ.tail_RegulationModes.csv")
BZNY.tail_RegulationModes <- read.csv("BZNY.tail_RegulationModes.csv")
NYBZ.liver_RegulationModes <- read.csv("NYBZ.liver_RegulationModes.csv")
BZNY.liver_RegulationModes <- read.csv("BZNY.liver_RegulationModes.csv")

mat_either_genes_tail <- read.csv("mat_either_genes_new.csv")
mat_either_genes_liver <- read.csv("mat_either_genes_liver_new.csv")

gen_genes_tail <- read.csv("gen_genes_tail.csv")
gen_genes_liver <- read.csv("gen_genes_liver.csv")



####TAIL analysis ####

#ALL GENES - tail
#all genes in regulatory data
allGenes_hybrid_merge_tail <- dplyr::full_join(BZNY.tail_RegulationModes,NYBZ.tail_RegulationModes, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC","parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(across(category_BZNY, replace_na, "conserved"), #change NAs to conserved
         across(category_NYBZ, replace_na, "conserved"))
  
#add a column that contains the regulatory type for both hybrid comparisons
allGenes_hybrid_merge_tail$both <- paste(allGenes_hybrid_merge_tail$category_BZNY, allGenes_hybrid_merge_tail$category_NYBZ, sep="_")

#summarize data in plottable format
allGenes_hybrid_merge_tail_formatted <- allGenes_hybrid_merge_tail %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_NYBZ = factor(category_NYBZ, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#define color palettes for the differend regulatory types
#plot of all genes alluvial between hybrids in tail
allGenes_tail_alluvial <- allGenes_hybrid_merge_tail_formatted %>%
  ggplot(aes(y = n, axis1 = category_NYBZ, axis2 = category_BZNY)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.99, width = 0.33, color = "transparent") + ylab(("Number of genes")) +
  theme_classic()  + scale_x_continuous(breaks = c(1, 2), labels = c("NYBZ"="SARAxMANA","BZNY"="MANAxSARA")) + theme(legend.position = "none",
                          axis.text.x = element_text(face="bold", size =9, colour = 'black'), 
                          axis.title.y = element_text(face="bold", size =11),
                          axis.text.y = element_text(size =8, colour = 'black'))
allGenes_tail_alluvial





#GENETIC GENES - tail
#merge genetic genes with regulatory data
BZNYregModes_gen_genes_tail <- na.omit(merge(BZNY.tail_RegulationModes,gen_genes_tail,by="gene", all.y = TRUE)) #265 have reg modes
NYBZregModes_gen_genes_tail <- na.omit(merge(NYBZ.tail_RegulationModes,gen_genes_tail,by="gene", all.y = TRUE)) #276 have reg modes

#merge
genGenes_hybrid_merge_tail <- dplyr::full_join(BZNYregModes_gen_genes_tail,NYBZregModes_gen_genes_tail, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC", "parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(across(category_BZNY, replace_na, "conserved"), #change NAs to conserved
         across(category_NYBZ, replace_na, "conserved"))

#add a column that contains the regulatory type for both hybrid comparisons
genGenes_hybrid_merge_tail$both <- paste(genGenes_hybrid_merge_tail$category_BZNY, genGenes_hybrid_merge_tail$category_NYBZ, sep="_")

#summarize data in plottable format
genGenes_hybrid_merge_tail_formatted <- genGenes_hybrid_merge_tail %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_BZNY = factor(category_BZNY, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#define color palettes for the differend regulatory types
#plot of all genes alluvial between hybrids in tail
genGenes_tail_alluvial <- genGenes_hybrid_merge_tail_formatted %>%
  ggplot(aes(y = n, axis1 = category_BZNY, axis2 = category_NYBZ)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.99, width = 0.33, color = "transparent") +
  theme_classic() + theme(legend.position = "none") + scale_x_continuous(breaks = c(1, 2), labels = c("BZNY", "NYBZ"))
genGenes_tail_alluvial



#MATERNAL GENES - tail
#merge materna genes with regulatory data
BZNYregModes_mat_genes_tail <- na.omit(merge(BZNY.tail_RegulationModes,mat_either_genes_tail,by="gene", all.y = TRUE)) #84 have reg modes, 53 in strict
NYBZregModes_mat_genes_tail <- na.omit(merge(NYBZ.tail_RegulationModes,mat_either_genes_tail,by="gene", all.y = TRUE)) #93 have reg modes, 56 in strict

#merge
matGenes_hybrid_merge_tail <- dplyr::full_join(BZNYregModes_mat_genes_tail,NYBZregModes_mat_genes_tail, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC", "parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(
    category_BZNY = tidyr::replace_na(category_BZNY, "conserved"),
    category_NYBZ = tidyr::replace_na(category_NYBZ, "conserved"))

#add a column that contains the regulatory type for both hybrid comparisons
matGenes_hybrid_merge_tail$both <- paste(matGenes_hybrid_merge_tail$category_BZNY, matGenes_hybrid_merge_tail$category_NYBZ, sep="_")

#summarize data in plottable format
matGenes_hybrid_merge_tail_formatted <- matGenes_hybrid_merge_tail %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_BZNY = factor(category_BZNY, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#define color palettes for the differend regulatory types
#plot of all genes alluvial between hybrids in tail
matGenes_tail_alluvial <- matGenes_hybrid_merge_tail_formatted %>%
  ggplot(aes(y = n, axis1 = category_BZNY, axis2 = category_NYBZ)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.99, width = 0.33, color = "transparent") +
  theme_classic() + theme(legend.position = "none") + scale_x_continuous(breaks = c(1, 2), labels = c("BZNY", "NYBZ"))
matGenes_tail_alluvial





####LIVER analysis ####

#ALL GENES - liver
#all genes in regulatory data
allGenes_hybrid_merge_liver <- dplyr::full_join(BZNY.liver_RegulationModes,NYBZ.liver_RegulationModes, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC", "parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(across(category_BZNY, replace_na, "conserved"), #change NAs to conserved
         across(category_NYBZ, replace_na, "conserved"))

#add a column that contains the regulatory type for both hybrid comparisons
allGenes_hybrid_merge_liver$both <- paste(allGenes_hybrid_merge_liver$category_BZNY, allGenes_hybrid_merge_liver$category_NYBZ, sep="_")

#summarize data in plottable format
allGenes_hybrid_merge_liver_formatted <- allGenes_hybrid_merge_liver %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_NYBZ = factor(category_NYBZ, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#define color palettes for the differend regulatory types
#plot of all genes alluvial between hybrids in liver
allGenes_liver_alluvial <- allGenes_hybrid_merge_liver_formatted %>%
  ggplot(aes(y = n, axis1 = category_NYBZ, axis2 = category_BZNY)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.95, width = 0.33, color = "transparent") + ylab(("Number of genes")) +
  theme_classic()  + scale_x_continuous(breaks = c(1, 2), labels = c("NYBZ"="SARAxMANA\n ","BZNY"="MANAxSARA\n ")) + theme(legend.position = "none",
                                                                                                                     axis.text.x = element_text(face="bold", size =9, colour = 'black'), 
                                                                                                                     axis.title.y = element_text(face="bold", size =11),
                                                                                                                     axis.text.y = element_text(size =8, colour = 'black'))
allGenes_liver_alluvial


#GENETIC GENES - liver
#merge genetic genes with regulatory data
BZNYregModes_gen_genes_liver <- na.omit(merge(BZNY.liver_RegulationModes,gen_genes_liver,by="gene", all.y = TRUE)) 
NYBZregModes_gen_genes_liver <- na.omit(merge(NYBZ.liver_RegulationModes,gen_genes_liver,by="gene", all.y = TRUE)) 

#merge
genGenes_hybrid_merge_liver <- dplyr::full_join(BZNYregModes_gen_genes_liver,NYBZregModes_gen_genes_liver, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC", "parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(across(category_BZNY, replace_na, "conserved"), #change NAs to conserved
         across(category_NYBZ, replace_na, "conserved"))

#add a column that contains the regulatory type for both hybrid comparisons
genGenes_hybrid_merge_liver$both <- paste(genGenes_hybrid_merge_liver$category_BZNY, genGenes_hybrid_merge_liver$category_NYBZ, sep="_")

#summarize data in plottable format
genGenes_hybrid_merge_liver_formatted <- genGenes_hybrid_merge_liver %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_BZNY = factor(category_BZNY, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#define color palettes for the differend regulatory types
#plot of all genes alluvial between hybrids in liver
genGenes_liver_alluvial <- genGenes_hybrid_merge_liver_formatted %>%
  ggplot(aes(y = n, axis1 = category_BZNY, axis2 = category_NYBZ)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.99, width = 0.33, color = "transparent") +
  theme_classic() + theme(legend.position = "none") + scale_x_continuous(breaks = c(1, 2), labels = c("BZNY", "NYBZ"))
genGenes_liver_alluvial




#MATERNAL GENES - liver
#merge materna genes with regulatory data
BZNYregModes_mat_genes_liver <- na.omit(merge(BZNY.liver_RegulationModes,mat_either_genes_liver,by="gene", all.y = TRUE)) #84 have reg modes, 53 in strict
NYBZregModes_mat_genes_liver <- na.omit(merge(NYBZ.liver_RegulationModes,mat_either_genes_liver,by="gene", all.y = TRUE)) #93 have reg modes, 56 in strict

#merge
matGenes_hybrid_merge_liver <- dplyr::full_join(BZNYregModes_mat_genes_liver,NYBZregModes_mat_genes_liver, by = "gene") %>%
  select("gene", "category_BZNY" = "regulation_mode.x", "category_NYBZ" = "regulation_mode.y", "BZNY.NYvBZ_Log2FC" ,"NYBZ.NYvBZ_Log2FC", "parent.NYvBZ_Log2FC.x") %>% #just get the columns we care about
  dplyr::mutate(
    category_BZNY = tidyr::replace_na(category_BZNY, "conserved"),
    category_NYBZ = tidyr::replace_na(category_NYBZ, "conserved"))

#add a column that contains the regulatory type for both hybrid comparisons
matGenes_hybrid_merge_liver$both <- paste(matGenes_hybrid_merge_liver$category_BZNY, matGenes_hybrid_merge_liver$category_NYBZ, sep="_")

#summarize data in plottable format
matGenes_hybrid_merge_liver_formatted <- matGenes_hybrid_merge_liver %>%
  dplyr::group_by(category_BZNY, category_NYBZ, both) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories - all the different possible combinations
  dplyr::mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  dplyr::filter(both != "conserved_conserved") %>%# removes the category that is just the genes that are conserved in both hybrids - which is most of them
  dplyr::mutate(category_BZNY = factor(category_BZNY, levels = c("conserved","compensatory","cis x trans","cis + trans", "trans only","cis only"))) #order them for plotting

#plot of all genes alluvial between hybrids in liver
matGenes_liver_alluvial <- matGenes_hybrid_merge_liver_formatted %>%
  ggplot(aes(y = n, axis1 = category_BZNY, axis2 = category_NYBZ)) +
  geom_alluvium(aes(fill= category_BZNY), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 0.3, aes.bind = "flows") +
  scale_fill_manual(breaks = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only"),
                    values = c("darkgray", "#3399ff", "#03c04a", "#ffae42","#9966cd","#f95738"),
                    labels = c("conserved", "compensatory", "cis x trans","cis + trans","trans only", "cis only")) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.99, width = 0.33, color = "transparent") +
  theme_classic() + theme(legend.position = "none") + scale_x_continuous(breaks = c(1, 2), labels = c("BZNY", "NYBZ"))
matGenes_liver_alluvial





####CIS AND TRANS ONLY fixed expression across hybrids####


#look at how genes in the regulatory types differ in their expression direction across hybrids (ie. higher in BZ allele vs NY) - FLAWED
#need to actually look at difference in sign between higher in BZNY whole mice vs NYBZ whole mice (differences between BZ and NY allele dont matter in trans genes)

#find BZNY vs. NYBZ differential expression
#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#read in data and metadata
#these files are all samples together, can later be broken up for specific comparisons
all_countData <- read.csv('all_samples_counts.csv', header = TRUE, sep = ",")
all_metadata <- read.csv('all_samples_metadata.csv', header = TRUE, sep = ",")
all_metadata$id <- gsub("-", ".", all_metadata$id)

#separate tissues
tail_countData <- all_countData[ , colnames(all_countData) %like% "tail" | colnames(all_countData) %like% "gene"]
tail_metadata <- subset(all_metadata, tissue == "tail")
liver_countData <- all_countData[ , colnames(all_countData) %like% "liver" | colnames(all_countData) %like% "gene"]
liver_metadata <- subset(all_metadata, tissue == "liver")

#retain just F1s
tailF1_countData <- tail_countData[ , colnames(tail_countData) %like% "NB" | colnames(tail_countData) %like% "BN" | colnames(tail_countData) %like% "gene"]
tailF1_metadata <- subset(tail_metadata, pop == "BZNY" | pop == "NYBZ")
liverF1_countData <- liver_countData[ , colnames(liver_countData) %like% "NB" | colnames(liver_countData) %like% "BN" | colnames(liver_countData) %like% "gene"]
liverF1_metadata <- subset(liver_metadata, pop == "BZNY" | pop == "NYBZ")

#run DESeq
ddsBZNYvNYBZ_tail <- DESeqDataSetFromMatrix(countData=tailF1_countData, 
                                  colData=tailF1_metadata, 
                                  design=~pop, tidy = TRUE)
keepBZNYvNYBZ_tail <- rowMeans(counts(ddsBZNYvNYBZ_tail)) >= 10
ddsBZNYvNYBZ_tail <- ddsBZNYvNYBZ_tail[keepBZNYvNYBZ_tail,] 
deseq_BZNYvNYBZ_tail <- DESeq(ddsBZNYvNYBZ_tail)
res_BZNYvNYBZ_tail <- results(deseq_BZNYvNYBZ_tail, contrast=c("pop","NYBZ","BZNY")) #neg lfc = higher in BZNY
BZNYvNYBZ_tail <- as.data.frame(res_BZNYvNYBZ_tail)

ddsBZNYvNYBZ_liver <- DESeqDataSetFromMatrix(countData=liverF1_countData, 
                                            colData=liverF1_metadata, 
                                            design=~pop, tidy = TRUE)
keepBZNYvNYBZ_liver <- rowMeans(counts(ddsBZNYvNYBZ_liver)) >= 10
ddsBZNYvNYBZ_liver <- ddsBZNYvNYBZ_liver[keepBZNYvNYBZ_liver,] 
deseq_BZNYvNYBZ_liver <- DESeq(ddsBZNYvNYBZ_liver)
res_BZNYvNYBZ_liver <- results(deseq_BZNYvNYBZ_liver, contrast=c("pop","NYBZ","BZNY")) #neg lfc = higher in BZNY
BZNYvNYBZ_liver <- as.data.frame(res_BZNYvNYBZ_liver)

#get sig genes
sigDEGsFC05_BZNYvNYBZ_tail <- BZNYvNYBZ_tail %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #1503, neg LFC = higher expression in BZNY 
sigDEGsFC05_BZNYvNYBZ_tail$gene <- row.names(sigDEGsFC05_BZNYvNYBZ_tail)
sigDEGsFC05_BZNYvNYBZ_liver <- BZNYvNYBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #1133, neg LFC = higher expression in BZNY 
sigDEGsFC05_BZNYvNYBZ_liver$gene <- row.names(sigDEGsFC05_BZNYvNYBZ_liver)



#TAIL

#merge significant DE dataframes with allGenes reg mode dataset to see how many are differentially expressed
allGenes_hybrid_merge_tail_sig <- na.omit(merge(allGenes_hybrid_merge_tail,sigDEGsFC05_BZNYvNYBZ_tail,by="gene", all.y = TRUE)) #130
#calculate number of sig diff shared cis across hybrids / overall cis shared across hybrids
(sum(allGenes_hybrid_merge_tail_sig$both == "cis only_cis only")/sum(allGenes_hybrid_merge_tail$both == "cis only_cis only")) #0.0156, 11/705
(sum(allGenes_hybrid_merge_tail_sig$both == "trans only_trans only")/sum(allGenes_hybrid_merge_tail$both == "trans only_trans only")) #0.04399, 41/932
#the robustness of cis and trans changes in the tail (all genes) - p-value = 0.001183 via chi sq test
prop.test(x = c(11, 41), n = c(705, 932), correct = FALSE)

#merge significant DE dataframes with matGenes reg mode dataset to see how many are differentially expressed
matGenes_hybrid_merge_tail_sig <- na.omit(merge(matGenes_hybrid_merge_tail,sigDEGsFC05_BZNYvNYBZ_tail,by="gene", all.y = TRUE)) #74
#calculate number of sig diff shared cis across hybrids / overmat cis shared across hybrids
(sum(matGenes_hybrid_merge_tail_sig$both == "cis only_cis only")/sum(matGenes_hybrid_merge_tail$both == "cis only_cis only")) #0.625, 5/8
(sum(matGenes_hybrid_merge_tail_sig$both == "trans only_trans only")/sum(matGenes_hybrid_merge_tail$both == "trans only_trans only")) #0.524, 32/61
#the robustness of cis and trans changes in the tail (mat genes) - p-value =0.716 via fishers test
fisher.test(matrix(c(5,3,32,29), nrow = 2))

#merge significant DE dataframes with genGenes reg mode dataset to see how many are differentially expressed
genGenes_hybrid_merge_tail_sig <- na.omit(merge(genGenes_hybrid_merge_tail,sigDEGsFC05_BZNYvNYBZ_tail,by="gene", all.y = TRUE)) #16
#calculate number of sig diff shared cis across hybrids / overall cis shared across hybrids
(sum(genGenes_hybrid_merge_tail_sig$both == "cis only_cis only")/sum(genGenes_hybrid_merge_tail$both == "cis only_cis only")) #0.028, 3/107
(sum(genGenes_hybrid_merge_tail_sig$both == "trans only_trans only")/sum(genGenes_hybrid_merge_tail$both == "trans only_trans only")) #0.032, 1/31
#the robustness of cis and trans changes in the tail (gen genes) - p-value = 1 via fishers test
fisher.test(matrix(c(3,104,1,30), nrow = 2))




#LIVER

#merge significant DE dataframes with allGenes reg mode dataset to see how many are differentially expressed
allGenes_hybrid_merge_liver_sig <- na.omit(merge(allGenes_hybrid_merge_liver,sigDEGsFC05_BZNYvNYBZ_liver,by="gene", all.y = TRUE)) #170
#calculate number of sig diff shared cis across hybrids / overall cis shared across hybrids
(sum(allGenes_hybrid_merge_liver_sig$both == "cis only_cis only")/sum(allGenes_hybrid_merge_liver$both == "cis only_cis only")) #0.04767064, 44/923
(sum(allGenes_hybrid_merge_liver_sig$both == "trans only_trans only")/sum(allGenes_hybrid_merge_liver$both == "trans only_trans only")) #0.0547619, 23/420
#the robustness of cis and trans changes in the liver (all genes) - p-value = 0.58 via chi sq test
prop.test(x = c(44, 23), n = c(923, 420), correct = FALSE)

#merge significant DE dataframes with matGenes reg mode dataset to see how many are differentially expressed
matGenes_hybrid_merge_liver_sig <- na.omit(merge(matGenes_hybrid_merge_liver,sigDEGsFC05_BZNYvNYBZ_liver,by="gene", all.y = TRUE)) #33
#calculate number of sig diff shared cis across hybrids / overmat cis shared across hybrids
(sum(matGenes_hybrid_merge_liver_sig$both == "cis only_cis only")/sum(matGenes_hybrid_merge_liver$both == "cis only_cis only")) #0.4, 8/20
(sum(matGenes_hybrid_merge_liver_sig$both == "trans only_trans only")/sum(matGenes_hybrid_merge_liver$both == "trans only_trans only")) #0.4667, 7/15
#the robustness of cis and trans changes in the liver (mat genes) - p-value = 0.7412 via fishers test
fisher.test(matrix(c(8,12,7,8), nrow = 2))

#merge significant DE dataframes with genGenes reg mode dataset to see how many are differentially expressed
genGenes_hybrid_merge_liver_sig <- na.omit(merge(genGenes_hybrid_merge_liver,sigDEGsFC05_BZNYvNYBZ_liver,by="gene", all.y = TRUE)) #57
#calculate number of sig diff shared cis across hybrids / overall cis shared across hybrids
(sum(genGenes_hybrid_merge_liver_sig$both == "cis only_cis only")/sum(genGenes_hybrid_merge_liver$both == "cis only_cis only")) #0.08, 22/275
(sum(genGenes_hybrid_merge_liver_sig$both == "trans only_trans only")/sum(genGenes_hybrid_merge_liver$both == "trans only_trans only")) #0.13, 11/82
#the robustness of cis and trans changes in the liver (gen genes) - p-value = 0.19 via fishers test
fisher.test(matrix(c(22,253,11,71), nrow=2))




#maternal trans genes restricted one hybrid - connected directly to trans regulated maternal effect
BZNYonly_matTrans_tail <- matGenes_hybrid_merge_tail %>% dplyr::filter(category_BZNY == "trans only" & category_NYBZ != "trans only") #23 genes
NYBZonly_matTrans_tail <- matGenes_hybrid_merge_tail %>% dplyr::filter(category_NYBZ == "trans only" & category_BZNY != "trans only") #12 genes

#copy these genes
library(clipr) #too small of gene lists to look at GO terms
write_clip(BZNYonly_matTrans_tail$gene) #in a simple GO there are no enrichments here
write_clip(NYBZonly_matTrans_tail$gene) #in a simple GO there are no enrichments here

#what about those trans genes that change expression between hybrids?
matGenes_hybrid_merge_tail_sig_tt <- matGenes_hybrid_merge_tail_sig %>% filter(matGenes_hybrid_merge_tail_sig$both == "trans only_trans only") #32
write_clip(matGenes_hybrid_merge_tail_sig_tt$gene) #no sig enrichment here





#### Plasticity cis-eQTL #####
#genes regulated in cis but only in one F1 - trans only is the other one
allGenes_cisBZNY_transNYBZ_liver <- allGenes_hybrid_merge_liver %>% filter(allGenes_hybrid_merge_liver$both == "cis only_trans only") #35 - no GO
allGenes_cisNYBZ_transBZNY_liver <- allGenes_hybrid_merge_liver %>% filter(allGenes_hybrid_merge_liver$both == "trans only_cis only") #31 - no GO

allGenes_cisBZNY_transNYBZ_tail <- allGenes_hybrid_merge_tail %>% filter(allGenes_hybrid_merge_tail$both == "cis only_trans only") #57 - no GO
allGenes_cisNYBZ_transBZNY_tail <- allGenes_hybrid_merge_tail %>% filter(allGenes_hybrid_merge_tail$both == "trans only_cis only") #60 - many exocytosis terms...



#hone in on some genes of interest

#liver - focus on BZNY cis because the BZ mother is inducing the large weight effect - it could be either actually
allGenes_cisBZNY_transNYBZ_liver_mat <- merge(allGenes_cisBZNY_transNYBZ_liver, mat_either_genes_liver, by = "gene")
allGenes_cisBZNY_transNYBZ_liver_parentDiff <- allGenes_cisBZNY_transNYBZ_liver %>% filter(gene %in% row.names(sigDEGsFC05_NNvBB_liver))#from LIVERmat_vs_gen.R
allGenes_cisBZNY_transNYBZ_liver_NYCFdiff <- allGenes_cisBZNY_transNYBZ_liver %>% filter(gene %in% row.names(sigDEGsFC05_NNvNNCF_liver))#from LIVERmat_vs_gen.R

allGenes_cisNYBZ_transBZNY_liver_mat <- merge(allGenes_cisNYBZ_transBZNY_liver, mat_either_genes_liver, by = "gene")
allGenes_cisNYBZ_transBZNY_liver_NYCFdiff <- allGenes_cisNYBZ_transBZNY_liver %>% filter(gene %in% row.names(sigDEGsFC05_NNvNNCF_liver))#from LIVERmat_vs_gen.R

NY_CFonlyGenes_regTypes <- merge(NY_CFonlyGenes, allGenes_hybrid_merge_liver, by="gene")
NY_CFonlyGenes_regTypes_MP <- merge(NY_CFonlyGenes_regTypes, BW_MPgenes_full, by="gene")














#looking into our maternal trans genes in the tail - I feel like these could be telling us something
#what is that group of genes that stays trans regulated but switches directions? consistenly controlled by maternal factors but they have opposite effects depending on who mom is

#trans mat changes in the tail that switch their directionality
transSignChange_matTail <- matGenes_hybrid_merge_tail_CT %>% filter(regSign == "trans only_NY high_trans only_BZ high" | regSign == "trans only_BZ high_trans only_NY high") #20 genes
transSignSame_matTail <- matGenes_hybrid_merge_tail_CT %>% filter(regSign == "trans only_NY high_trans only_NY high" | regSign == "trans only_BZ high_trans only_BZ high") #41 genes
transOnly_matTail <- matGenes_hybrid_merge_tail_CT %>% filter(both == "trans only_trans only") #61 genes

#copy these 20 genes
library(clipr)
write_clip(transSignChange_matTail$gene) #in a simple GO there are no enrichments here

#are these sign change trans genes like...imprinted??

#ENSMUSG00000001494 - sost, there are bone phenotypes
#ENSMUSG00000005360 - Slc1a3, central nervous system
#ENSMUSG00000006235 - Epor, erythropoietin - production of red blood cells
#ENSMUSG00000008843 - Cldn13 - mouse specific, related to red blood cells
#ENSMUSG00000019194 - Scn1b - has growth effects, neuronal effects
#ENSMUSG00000024987 - Cyp26a1 - limb effects
#ENSMUSG00000025889 - Snca, coordination
#ENSMUSG00000026885 - Ttll11, hematopoeitic
#ENSMUSG00000026896 - Ifih1, hematopoietic, immune
#ENSMUSG00000027510 - Rbm38, hematopoietic
#ENSMUSG00000027778 - Ift80, short bones when knocked out
#ENSMUSG00000027887 - Sypl2, growth and skeletal muscle effects
#ENSMUSG00000028332 - Hemgn, also related to red blood cells and growth
#ENSMUSG00000029862 - Clcn1, spinal abnormalities and body size variation
#ENSMUSG00000030592 - Ryr1, skeletal development effects
#ENSMUSG00000032131 - Abcg4, lipids and cholesterol
#ENSMUSG00000036854 - Hspb6, a heat shock protein...if BZ mom it is high in BZ and if NY mom it is high in NY
#ENSMUSG00000043065 - Spice1
#ENSMUSG00000043795
#ENSMUSG00000051839


write_clip(transOnly_matTail$gene) #in a simple GO there are no interesting enrichments here
write_clip(transSignSame_matTail$gene) #in a simple GO there are no interesting enrichments here
