#genetic vs. maternal control for NY vs BZ liver gene expression differences

#defining maternally controlled genes and genetically controlled genes
#males and females together
#no F1s in this analysis
#only liver data

#started 1.26.26
#using results data with correct design: made in LIVERdefine_DE_comps.R

# 2.10.26 update: removing the third filter level of maternal genes because it removed mom specific expression changes 
#(required a sig. maternal state DESeq test which reflects genes where both NY mom and BZ mom are driving a difference)
#also focusing on the 0.5 lfc set of genes

#2.12.26 update: introducing a strict maternal gene definition which requires there NOT to be a difference between pups with the same mom (ie BZ-CF vs NY for NY dam and NY-CF vs BZ for BZ dam)

#load packages
library(DESeq2)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggrepel)

#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_DE_comps.RData")

#Define significant DEGs, going to add in log2FC of 0.5 from the start because there are so many diff exp genes at just p<0.05 (10,000 ish)
sigDEGsFC05_NNvBB_liver <- pureline_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3871, neg LFC = higher expression in BZ 
sigDEGsFC05_NNCFvBBCF_liver <- cf_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #4595, neg LFC = higher expression in BZ 
sigDEGsFC05_NNvNNCF_liver <- NY_pureVcf_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #836, neg LFC = higher expression NY-CF
sigDEGsFC05_BBCFvBB_liver <- BZ_cfVpure_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #309, neg LFC = higher expression BZ-PURELINE!!
sigDEGsFC05_NNCFvBB_liver <- NYCFvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #4329, neg LFC = higher expression in BZ pureline
sigDEGsFC05_NNvBBCF_liver <- NYvsBZCF_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3810, neg LFC = higher expression BZ-CF!!
sigDEGsFC05_genoG_liver <- liver_G %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #3982, neg LFC = higher expression in BZ
sigDEGsFC05_maternalE_liver <- liver_E %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #465, neg LFC = higher expression in BZ reared mice
sigDEGsFC05_GxE_liver<- liver_GxE %>% filter(padj<=0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) #427

#Define significant DEGs, going to add in log2FC of 0.25 from the start because there are so many diff exp genes at just p<0.05, but 0.5 is a bit much
sigDEGsFC25_NNvBB_liver <- pureline_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_NNCFvBBCF_liver <- cf_NYvsBZ_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_NNvNNCF_liver <- NY_pureVcf_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_BBCFvBB_liver <- BZ_cfVpure_liver %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_genoG_liver <- liver_G %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_maternalE_liver <- liver_E %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 
sigDEGsFC25_GxE<- liver_GxE %>% filter(padj<=0.05 & (log2FoldChange > 0.25 | log2FoldChange < -0.25)) 



#### Maternal sets of genes ####

#For Brazil dams, a maternal set of genes involved in liver expression variation is: 
#different between NYNY v BZBZ, 
#different between NYNY v NYNY-CF (because NYNY-CF is reared by BZ)
#and has a significant environmental component through the "E" DESeq comparison (maternal status significantly influenced expression)  - REMOVED in 2.10.26 update
#mat_BZ_genes_liver <- sigDEGsFC05_NNvBB_liver %>% filter(row.names(.) %in% row.names(sigDEGsFC05_NNvNNCF_liver) & row.names(.) %in% row.names(sigDEGsFC05_maternalE_liver)) #217 genes

#different between NYNY v BZBZ, 
#different between NYNY v NYNY-CF (because NYNY-CF is reared by BZ)
#sign is the same in the two comparisons - mom is driving expression in the same direction
mat_BZ_genes_liver <- sigDEGsFC05_NNvBB_liver %>% rownames_to_column("gene") %>% 
  inner_join(sigDEGsFC05_NNvNNCF_liver %>% rownames_to_column("gene"), by = "gene") %>% filter(sign(log2FoldChange.x) == sign(log2FoldChange.y)) #257 genes (vs. 217 before)

#strict maternal filter - no difference between BZBZ and NYNY-CF
mat_BZ_genes_liver_strict <- subset(mat_BZ_genes_liver, !(mat_BZ_genes_liver$gene)%in%(row.names(sigDEGsFC05_NNCFvBB_liver))) #144 down from 257


#For New York dams, a maternal set of genes involved in liver length variation is: 
#different between NYNY v BZBZ, 
#different between BZBZ v BZBZ-CF (because BZBZ-CF is reared by NY)
#and has a significant environmental component through the "E" DESeq comparison (maternal status significantly influenced expression)  - REMOVED in 2.10.26 update
#mat_NY_genes_liver <- sigDEGsFC05_NNvBB_liver %>% filter(row.names(.) %in% row.names(sigDEGsFC05_BBCFvBB_liver)) & row.names(.) %in% row.names(sigDEGsFC05_maternalE_liver)) #111 genes - 

#different between NYNY v BZBZ, 
#different between NYNY v NYNY-CF (because NYNY-CF is reared by BZ)
#sign is the same in the two comparisons - mom is driving expression in the same direction
mat_NY_genes_liver <- sigDEGsFC05_NNvBB_liver %>% rownames_to_column("gene") %>% 
  inner_join(sigDEGsFC05_BBCFvBB_liver %>% rownames_to_column("gene"), by = "gene") %>% filter(sign(log2FoldChange.x) == sign(log2FoldChange.y)) #161 genes (vs. 111 before)

#strict maternal filter - no difference between BZBZ and NYNY-CF
mat_NY_genes_liver_strict <- subset(mat_NY_genes_liver, !(mat_NY_genes_liver$gene)%in%(row.names(sigDEGsFC05_NNvBBCF_liver))) #85 down from 161


#but do the mat lists overlap?
length(intersect(as.list(rownames(mat_BZ_genes_liver)), as.list(rownames(mat_NY_genes_liver)))) #161, 74 before - sorta
shared_mat_genes_liver <- merge(mat_BZ_genes_liver, mat_NY_genes_liver, by=0) #order: NNvBB, NNvNNCF, NNvBB, BBCFvBB

#merge BZ and NY maternal set to have a set of maternal in either genes
mat_either_genes_liver <- bind_rows(mat_BZ_genes_liver, mat_NY_genes_liver) %>% distinct(gene, .keep_all = TRUE) #355 vs 254 before
write.csv(mat_either_genes_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/mat_either_genes_liver_new.csv")

mat_either_genes_liver_strict <- bind_rows(mat_BZ_genes_liver_strict, mat_NY_genes_liver_strict) %>% distinct(gene, .keep_all = TRUE) #206 vs 355 before
write.csv(mat_either_genes_liver_strict, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/mat_either_genes_liver_strict.csv")



#### Genetic sets of genes ####

#Just one set among NY and BZ
#different between NYNY v BZBZ, 
#different BZBZ-CF v NYNY-CF, 
#and has a significant genotype component through the "G" DESeq comparison (genotype status significantly influenced expression)
#also filter for genes where the direction is the same between NYvBZ pureline and NYvBZ CF
gen_genes_liver <- merge(sigDEGsFC05_NNvBB_liver, sigDEGsFC05_NNCFvBBCF_liver, by=0) #2866 #many more than tail (1000 more)
gen_genes_liver <- subset(gen_genes_liver, Row.names %in% row.names(sigDEGsFC05_genoG_liver)) #2812
gen_genes_liver_misExp <- gen_genes_liver %>% dplyr::filter(sign(log2FoldChange.x) != sign(log2FoldChange.y)) #2798 - 14 are misexpressed
gen_genes_liver <- gen_genes_liver %>% filter(sign(gen_genes_liver$log2FoldChange.x) == sign(gen_genes_liver$log2FoldChange.y)) #2798 - 14 are misexpressed
colnames(gen_genes_liver)[1] <- "gene"

write.csv(gen_genes_liver, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/gen_genes_liver.csv")



#### Venn Diagrams ####

#pureline vs CF NYvBZ differences  -  significant genes (removing those that are not in the same direction for the overlap)
#remove the full overlap from the comparison, but have the the overlap venn size just be those that are in the same direction
library(eulerr)

plot(euler(c("NY v BZ \npureline"=(3871-2812),"NY v BZ \ncross\nfostered"=(4595-2812),
             "NY v BZ \npureline&NY v BZ \ncross\nfostered"=2812)), quantities = TRUE,
     fills=c("#00BFC4","#F8766D","#C77CFF"), alpha = 0.7)




#pureline NYvBZ and maternal set overlaps

#NY dams
length(intersect(row.names(sigDEGsFC05_NNvBB_liver), (mat_NY_genes_liver$gene))) #161, was 111
plot(euler(c("NY v BZ \npureline"=(3871-161),"NY maternal \ngenes"=(0),
             "NY v BZ \npureline&NY maternal \ngenes"=161)), quantities = TRUE,
     fills=c("#00BFC4","#F8766D"), alpha = 0.7)
#BZ dams
length(intersect(row.names(sigDEGsFC05_NNvBB_liver), (mat_BZ_genes_liver$gene))) #257, was 217
plot(euler(c("NY v BZ \npureline"=(3871-257),"BZ maternal \ngenes"=(0),
             "NY v BZ \npureline&BZ maternal \ngenes"=257)), quantities = TRUE,
     fills=c("#ce8c14","#F8766D"), alpha = 0.7)



##### Genes related to the maternal effect on liver ONLY in NY #####

#Genes with a unidirectional maternal effect in NYCF - so maternal genes induced by BZ mom, and then not different between BZ and BZCF
NY_CFonlyGenes <- mat_BZ_genes_liver %>% filter(!gene %in% row.names(sigDEGsFC05_BBCFvBB_liver))
NY_CFonlyGenes_MPgenes <- merge(NY_CFonlyGenes, BW_MPgenes_full, by="gene")



#specific genes 
#GxE examples - CF misexpression - 

#Aldoc - weird this just seems maternal - yeah it is in the maternal set
ggplot(plotCounts(dds_liver, gene = "ENSMUSG00000017390", intgroup = "pop", returnData=TRUE), aes(x=pop, y=count, color=pop)) +
  geom_point(position=position_jitter(w=0.1,h=0), size=2) + ggtitle("Gene Expression of Aldoc") + theme_bw() #related to cholesterol

#aquaporin 3 - high in NY-C only - starvation related? this is in our BZ maternal genes
ggplot(plotCounts(dds_liver, gene = "ENSMUSG00000028435", intgroup = "pop", returnData=TRUE), aes(x=pop, y=count, color=pop)) +
  geom_point(position=position_jitter(w=0.1,h=0), size=2) + ggtitle("Gene Expression of aquaporin 3") + theme_bw() #related to water intake

#	apolipoprotein A-IV - high in NY-C only - starvation related? this is in our BZ maternal genes
ggplot(plotCounts(dds_liver, gene = "ENSMUSG00000032080", intgroup = "pop", returnData=TRUE), aes(x=pop, y=count, color=pop)) +
  geom_point(position=position_jitter(w=0.1,h=0), size=2) + ggtitle("Gene Expression of 	apolipoprotein A-IV") + theme_bw() #related to cholesterol
