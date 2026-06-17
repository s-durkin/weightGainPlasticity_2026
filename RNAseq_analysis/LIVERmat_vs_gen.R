#code to determine genetic vs. maternal control for NY vs BZ liver gene expression differences

#Corresponds to patterns shown in Figure 2D

#defining maternally controlled genes and genetically controlled genes
#males and females together
#no F1s in this analysis
#only liver data


#load packages
library(DESeq2)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggrepel)

#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_DE_comps.RData")

#Define significant DEGs, p<0.05 and log2FC of +- 0.5 
sigDEGsFC05_NNvBB_liver <- pureline_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3871, neg LFC = higher expression in BZ 
sigDEGsFC05_NNCFvBBCF_liver <- cf_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #4595, neg LFC = higher expression in BZ 
sigDEGsFC05_NNvNNCF_liver <- NY_pureVcf_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #836, neg LFC = higher expression NY-CF
sigDEGsFC05_BBCFvBB_liver <- BZ_cfVpure_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #309, neg LFC = higher expression BZ-PURELINE!!
sigDEGsFC05_NNCFvBB_liver <- NYCFvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #4329, neg LFC = higher expression in BZ pureline
sigDEGsFC05_NNvBBCF_liver <- NYvsBZCF_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3810, neg LFC = higher expression BZ-CF!!
sigDEGsFC05_genoG_liver <- liver_G %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3982, neg LFC = higher expression in BZ
sigDEGsFC05_maternalE_liver <- liver_E %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #465, neg LFC = higher expression in BZ reared mice
sigDEGsFC05_GxE_liver<- liver_GxE %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #427





#### Maternal sets of genes ####

#For Brazil dams, a maternal set of genes involved in liver expression variation is: 
#different between NYNY v BZBZ, 
#different between NYNY v NYNY-CF (because NYNY-CF is reared by BZ)
#sign is the same in the two comparisons - mom is driving expression in the same direction
mat_BZ_genes_liver <- sigDEGsFC05_NNvBB_liver %>% rownames_to_column("gene") %>% 
  inner_join(sigDEGsFC05_NNvNNCF_liver %>% rownames_to_column("gene"), by = "gene") %>% filter(sign(log2FoldChange.x) == sign(log2FoldChange.y)) #257 genes 


#For New York dams, a maternal set of genes involved in liver length variation is: 
#different between NYNY v BZBZ, 
#different between NYNY v NYNY-CF (because NYNY-CF is reared by BZ)
#sign is the same in the two comparisons - mom is driving expression in the same direction
mat_NY_genes_liver <- sigDEGsFC05_NNvBB_liver %>% rownames_to_column("gene") %>% 
  inner_join(sigDEGsFC05_BBCFvBB_liver %>% rownames_to_column("gene"), by = "gene") %>% filter(sign(log2FoldChange.x) == sign(log2FoldChange.y)) #161 genes 


#merge BZ and NY maternal set to have a set of maternal in either genes
mat_either_genes_liver <- bind_rows(mat_BZ_genes_liver, mat_NY_genes_liver) %>% distinct(gene, .keep_all = TRUE) #355 genes
write.csv(mat_either_genes_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/mat_either_genes_liver_new.csv")





#### Genetic sets of genes ####

#Just one set among NY and BZ
#different between NYNY v BZBZ, 
#different BZBZ-CF v NYNY-CF, 
#and has a significant genotype component through the "G" DESeq comparison (genotype status significantly influenced expression)
#also filter for genes where the direction is the same between NYvBZ pureline and NYvBZ CF
gen_genes_liver <- merge(sigDEGsFC05_NNvBB_liver, sigDEGsFC05_NNCFvBBCF_liver, by=0) #2866 
gen_genes_liver <- subset(gen_genes_liver, Row.names %in% row.names(sigDEGsFC05_genoG_liver)) #2812
gen_genes_liver_misExp <- gen_genes_liver %>% dplyr::filter(sign(log2FoldChange.x) != sign(log2FoldChange.y)) #2798 - 14 are misexpressed
gen_genes_liver <- gen_genes_liver %>% filter(sign(gen_genes_liver$log2FoldChange.x) == sign(gen_genes_liver$log2FoldChange.y)) #2798 - 14 are misexpressed
colnames(gen_genes_liver)[1] <- "gene"

write.csv(gen_genes_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/gen_genes_liver.csv")



