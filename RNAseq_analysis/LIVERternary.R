#liver ternary plots for displaying G, E, and GxE (Figure 1F)

#referring to Ballinger et al. 2023 PNAS github: https://github.com/malballinger/BallingerMack_PNAS_2023/blob/main/code/postprocess_RNAseq/plot_parents_ternary_GxE.R

#uses datasets defined in LIVERdefine_DE&ASE_comps.R


#### set-up ####

#load packages
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_DE_comps.RData")

#GxE lists

# collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_GxE_liver <- liver_GxE %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
# save a list of GxE genes defined by DESEQ
DESEQ_GxE_list_liver <- DESEQ_GxE_liver %>% pull(gene) #427



#### prep data ####


#rename column in G, E and GxE data
base_resG_liver <- liver_G %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(abslfc = abs(log2FoldChange))

base_resE_liver <- liver_E %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(abslfc = abs(log2FoldChange))

base_resGxE_liver <- liver_GxE %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(abslfc = abs(log2FoldChange))

merge_resGE_liver <- dplyr::full_join(base_resG_liver, base_resE_liver, by = "gene") %>% 
  dplyr::select("gene","lfc.g" = "log2FoldChange.x","padj.g" = "padj.x","abslfc.g" = "abslfc.x","lfc.e" = "log2FoldChange.y","padj.e" = "padj.y","abslfc.e" = "abslfc.y")

mergeall_resGE_liver <- dplyr::full_join(merge_resGE_liver, base_resGxE_liver, by = "gene") %>%
  dplyr::select("gene","lfc.g", "padj.g", "abslfc.g","lfc.e", "padj.e", "abslfc.e","lfc.gxe" = "log2FoldChange","padj.gxe" = "padj","abslfc.gxe" = "abslfc")


sig_groups_liver <- mergeall_resGE_liver %>%
  filter((padj.g < 0.05 & abslfc.g > 0.5) | (padj.e < 0.05 & abslfc.e >0.5) | (abslfc.gxe > 0.5) & (gene %in% DESEQ_GxE_list_liver)) %>%
  mutate(main_category = case_when(
    padj.g < 0.05 & (!gene %in% DESEQ_GxE_list_liver) ~ 1, # g_liver = 1
    padj.e < 0.05 & (!gene %in% DESEQ_GxE_list_liver) ~ 2, # e_liver = 2
    gene %in% DESEQ_GxE_list_liver ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
#4296 genes





#### ternary plot ####

# plot parental expression patterns of parents as a ternary plot
plot_liver <- ggtern(data = sig_groups_liver, aes(x=g_propvar, y=gxe_propvar, z=e_propvar)) +
  geom_point(fill = "black", color = "black", size = 1, stroke = 1, shape = 21, alpha = 0.12) +
  scale_L_continuous(breaks=c(.25,.50,.75,1)) +
  scale_R_continuous(breaks=c(.25,.50,.75,1)) +
  scale_T_continuous(breaks=c(.25,.50,.75,1)) +
  theme_bw(base_size=6)+
  theme_showarrows() +
  theme(axis.text = element_text(size = 10, family = "sans"),
        axis.title = element_text(size = 12, family = "sans"),
        aspect.ratio = 1/1,
        tern.axis.arrow = element_line(color = "black", size = 1)) +
  labs(x = "G",
       xarrow = NULL,
       y = "GxE",
       yarrow = NULL,
       z = "E",
       zarrow = NULL)
plot_liver







#### effects ####

n_total_genes_liver <- nrow(mergeall_resGE_liver) #16659
n_sig_genes_liver <- mergeall_resGE_liver %>% filter(padj.g < 0.05 & abslfc.g > 0.5 | padj.e < 0.05 & abslfc.e > 0.5| padj.gxe < 0.05 & abslfc.gxe > 0.5) %>% nrow()
#4296 with log fold change filter at 0.5

# effect sizes
# How many sig. DEG for genotype only, with no environmental influence?
g_liver <- mergeall_resGE_liver %>%
  filter((padj.g < 0.05 & abslfc.g > 0.5) &
           (!gene %in% DESEQ_GxE_list_liver)) # this list is generated above
n_g_genes_liver <- nrow(g_liver) #3747 at lfc filter of 0.5
effectsize_g_liver <- g_liver %>% summarise(effectsize = mean(abslfc.g)) #1.34 when log fold change filter is 0.5 

# How many sig. DEG for env only, with no genotype influence?
e_liver <- mergeall_resGE_liver %>%
  filter((padj.e < 0.05 & abslfc.e > 0.5) &
           (!gene %in% DESEQ_GxE_list_liver))
n_e_genes_liver <- nrow(e_liver) #346 when log fold change filter is 0.5 
effectsize_e_liver <- e_liver %>% summarise(effectsize = mean(abslfc.e)) #0.76 when log fold change filter is 0.5



#effects over all genes, not just those differentially G or E
effectsizeAll_g_liver <- sig_groups_liver %>% summarise(effectsize = mean(abslfc.g)) #1.25 when log fold change filter is 0.5
effectsizeAll_e_liver <- sig_groups_liver %>% summarise(effectsize = mean(abslfc.e)) #0.229 when log fold change filter is 0.5 



