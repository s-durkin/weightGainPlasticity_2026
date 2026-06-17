#Code to compare regulatory categories defined in NYBZ vs BZNY F1 hybrids. Figure 2F.


#Load packages
library(dplyr)
library(ggplot2)
library(tidyverse)

#read in Regulatory data for NYBZ and BZNY and the genetic and maternal genes
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")

#load regulatory modes spreadsheets - both tissues - generated in LIVERdefine_cis_trans.R
NYBZ.tail_RegulationModes <- read.csv("NYBZ.tail_RegulationModes.csv")
BZNY.tail_RegulationModes <- read.csv("BZNY.tail_RegulationModes.csv")
NYBZ.liver_RegulationModes <- read.csv("NYBZ.liver_RegulationModes.csv")
BZNY.liver_RegulationModes <- read.csv("BZNY.liver_RegulationModes.csv")

mat_either_genes_liver <- read.csv("mat_either_genes_liver_new.csv")

gen_genes_liver <- read.csv("gen_genes_liver.csv")




####Regulatory category stability analysis ####

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


