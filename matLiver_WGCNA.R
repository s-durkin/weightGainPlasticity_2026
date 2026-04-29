#WGCNA for  liver

#gene lists and dataframes from LIVERmat_vs_gen.R

#### load packages ####
library(WGCNA)
library(DESeq2)


#Navigate to correct folder. 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_DE_comps.RData")



#normalize counts from dds and convert to matrix
dds_liver_norm <- vst(dds_liver, blind = FALSE)
normalized_dds_liver <- assay(dds_liver_norm) 

#keep only those that are expressed and variable - WGCNA likes that
keep <- rowSums(normalized_dds_liver > 0) >= (0.5 * ncol(normalized_dds_liver))
normalized_dds_liver <- normalized_dds_liver[keep, ]
var_genes <- apply(normalized_dds_liver, 1, var)
normalized_dds_liver <- normalized_dds_liver[order(var_genes, decreasing = TRUE)[1:8000], ]

# transform the data - idk why but this doesnt work, but gets fized at mat gene filtering step
normalized_dds_liver <- t(normalized_dds_liver) 

#checking for good samples and genes - all good
gsg_liver <- goodSamplesGenes(normalized_dds_liver, verbose = 3)
normalized_dds_liver <- normalized_dds_liver[gsg_liver$goodSamples, gsg_liver$goodGenes]



#MATERNAL GENES - split before normalization
mat_normalized_dds_liver <- normalized_dds_liver[,colnames(normalized_dds_liver) %in% (as.list(mat_BZ_genes_liver$gene))] #from LIVERmat_vs_gen.R

#determine WGCNA parameters - mainly the power threshold
sft_matLiver <- pickSoftThreshold(mat_normalized_dds_liver, dataIsExpr = TRUE, powerVector = c(1:20), verbose = 5)

#calculate a measure of model fit for the threshold objects - R^2
sft_matLiver_df <- data.frame(sft_matLiver$fitIndices) %>% dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

#maternal genes: POWER = 8
ggplot(sft_matLiver_df, aes(x = Power, y = model_fit, label = Power)) + geom_point() + geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") + ylim(c(min(sft_matLiver_df$model_fit), 1.05)) + xlab("Soft Threshold (power)") + 
  ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale independence") + theme_classic()

matLiver_net <- blockwiseModules(mat_normalized_dds_liver,
                                TOMType = "signed", # topological overlap matrix, signed is better for RNAseq so modules reflect coactivation
                                power = 8, # soft threshold for network construction
                                minModuleSize = 20, #for gene sets 300-500 use 20, 500-2k use 30
                                randomSeed = 1234,# there's some randomness associated with this calculation
                                deepSplit = 2,
                                mergeCutHeight = 0.25
)

moduleColors_matLiver <- matLiver_net$colors
#check module size
table(moduleColors_matLiver) #blue has 66, grey 7, and turquoise 169 - blue is most correlated with maternal status

mod_genes_mat <- moduleColors_matLiver == "blue"
cor_mat_mat <- cor(mat_normalized_dds_liver[, mod_genes_mat])
mean(cor_mat_mat[upper.tri(cor_mat_mat)]) #mean corr is low: 0.01 - very low but its correlation with maternal status is so high (0.88) that it is still the best

mod_genes_mat <- moduleColors_matLiver == "grey"
cor_mat_mat <- cor(mat_normalized_dds_liver[, mod_genes_mat])
mean(cor_mat_mat[upper.tri(cor_mat_mat)]) #mean corr is 0.159 higher, but not really

mod_genes_mat <- moduleColors_matLiver == "turquoise"
cor_mat_mat <- cor(mat_normalized_dds_liver[, mod_genes_mat])
mean(cor_mat_mat[upper.tri(cor_mat_mat)]) #mean corr is low: 0.006


#going to try and relate to the trait of maternal group - I think WGCNA is really only gonna be interesting for our maternal liver genes here which is fine
#define a trait of maternal rearing group
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")
all_metadata <- read.csv('all_samples_metadata.csv', header = TRUE, sep = ",")
all_metadata$id <- gsub("-", ".", all_metadata$id)
liver_metadata <- subset(all_metadata, tissue == "liver" & pop != "BZNY" & pop != "NYBZ")

#convert trait of interest to numeric to run tests
mat_trait_liver <- model.matrix(~ 0 + maternal, data = liver_metadata)

#correlate maternal type with our module eigengenes
matLiver_MEs <- matLiver_net$MEs
moduleTraitCor_liver <- cor(matLiver_MEs, mat_trait_liver, use = "p")
moduleTraitP_liver <- corPvalueStudent(moduleTraitCor_liver, nSamples = nrow(mat_normalized_dds_liver))

#BLUE has strongest correlation with maternal types - lets zoom in - 0.88 corr with maternal type
moduleGenes_highSig_mat_liver <- moduleColors_matLiver == "blue"

#HUB GENES
#1st criteria: Module membership (kME) - How strongly a gene correlates with its module eigengene
kME_liver <- cor(mat_normalized_dds_liver[, moduleGenes_highSig_mat_liver], matLiver_MEs[, paste0("ME", "blue")], use = "p")

#2nd criteria: Gene significance (GS) - Gene Correlation with a trait (here mat group)
GS_liver <- cor(mat_normalized_dds_liver[, moduleGenes_highSig_mat_liver], mat_trait_liver[,"maternalBZ"], use = "p")

#high GS and high kME are often hub candidates
kME_liver <- as.vector(kME_liver)
GS_liver  <- as.vector(GS_liver)

names(kME_liver) <- colnames(mat_normalized_dds_liver)[moduleGenes_highSig_mat_liver]
names(GS_liver)  <- colnames(mat_normalized_dds_liver)[moduleGenes_highSig_mat_liver]

#plotting blue module relation to maternal status
boxplot(matLiver_MEs$MEblue ~ liver_metadata$maternal) #there is clear separation
hist(abs(GS_liver)) #most centered around 0.75 - high gene-trait correlation even though we have lower gene-gene correlation

#get a list of all blue module genes
blueModule_genes <- names(kME_liver)



#### HUB GENES ####

#defining hub genes
#want genes with high individual correlation with trait (GS) and correlation with module (kME)
hub_score <- abs(kME_liver) * abs(GS_liver)
top_hubs <- names(sort(hub_score, decreasing = TRUE))[1:10]
#Cyp2c55, Cyp3a11, Aass, Gpt, Sel1l3,Gm26793, Cth, Hmgcs1 (cholesterol),Suox (cysteine), Mvd (mevalonate)
#cysteine is implicated in catch-up growth in mice

#pick out high ones manually
hubCandidates_matLiver <- names(which(abs(kME_liver) > 0.8 & abs(GS_liver) > 0.4)) #29 - GO indicates many metabolic processes, specifically chloesterol (via precursors, saved as liverMatHubGO.csv)
hubCandidates_matLiver_strict <- names(which(abs(kME_liver) > 0.9 & abs(GS_liver) > 0.6)) #7
#pcsk9 (cholesterol), Gpt (glucose metabolism), Cyp2c55, Aass (lysine degredation), Sel1l3, Cyp3a11, Suox
#(there are two cytochrome P450 genes - genes in this complex are often linked to metabolism and weight differences)

write.csv(hubCandidates_matLiver,"~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/WGCNA_hubCandidates_matLiver.csv")
write.csv(hubCandidates_matLiver_strict,"~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/WGCNA_hubCandidates_matLiver_strict.csv")

#hub gene overlap with MP codes
intersect(hubCandidates_matLiver, as.list(BW_MPgenes_full$gene)) #from BW_QTLxLiverCis.R
#Il1rn, Cth, Aqp3(glycerol transport) (Il1rn and Aqp3 high in NY-CF, Cth strongly maternal)

#overlap with NY_CF only genes
intersect(hubCandidates_matLiver, as.list(NY_CFonlyGenes$gene)) #from NYBZvBZNY_regulatoryTypes.R, 10 of these
#Rgs16 (has diet implications), Ptafr, Il1rn, Rdh16 (depletion leads to fat accumulation, but NYCF very low), Aqp3
#Noxa1 (NADPH activator), D630039A03Rik, Cyp3a11, Bcas1, Suox, gm13387

#get regular names
hubGenes_named <- mm10_gene_clean %>% filter(gene %in% hubCandidates_matLiver) #from BW_QTLxLiverCis.R
modGenes_named <- mm10_gene_clean %>% filter(gene %in% blueModule_genes)



#### focus on just the genes related to BZ mothers because they are inducing the maternal effect on weight in NY-CF pups
#use pop not maternal type as the correlating variable and try to find modules where NY-CF is falling out

#convert trait of interest to numeric to run tests
pop_trait_liver <- model.matrix(~ 0 + pop, data = liver_metadata)

#correlate maternal type with our module eigengenes
popLiver_MEs <- matLiver_net$MEs
moduleTraitCor_popliver <- cor(popLiver_MEs, pop_trait_liver, use = "p")
moduleTraitP_popliver <- corPvalueStudent(moduleTraitCor_popliver, nSamples = nrow(mat_normalized_dds_liver))

#this also points to blue being most correlated with NY-CF





#### PLOTTING ####

#plot hub genes of interest
library(ggplot2)
library(dplyr)
library(purrr)

# Combine plotCounts for all genes
plot_data <- map_df(top_hubs, function(gene_id) {
  plotCounts(dds_liver,gene = gene_id,intgroup = "pop",returnData = TRUE) %>% mutate(gene = gene_id)})
plot_data$gene <- factor(plot_data$gene, levels = top_hubs)


#actual plot
custom_labels <- c("ENSMUSG00000025002" = "Cyp2c55","ENSMUSG00000056035" = "Cyp3a11","ENSMUSG00000029695" = "Aass",
                   "ENSMUSG00000022546" = "Gpt","ENSMUSG00000029189" = "Sel1l3","ENSMUSG00000032313" = "Gm26793",
                   "ENSMUSG00000028179" = "Cth","ENSMUSG00000093930" = "Hmgcs1","ENSMUSG00000049858" = "Suox","ENSMUSG00000006517" = "Mvd")

paletteGroup <- c("#ce8c14","#f2cd8b", "#006d92", "#7fbcd1")

ggplot(plot_data, aes(x = pop, y = count, color = pop)) + scale_color_manual(values=paletteGroup)+
  geom_jitter(width = 0.1, size = 1) +geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y", labeller = labeller(gene = custom_labels)) + labs(y = "Normalized counts", x = NULL) + 
  scale_x_discrete(labels = c("BZBZ" = "BZ","BZBZ-C" = "BZ-CF","NYNY" = "NY", "NYNY-C" = "NY-CF")) +
  theme_bw() + theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(size = 8), 
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")) 



#GO breakdown of the 29 hub genes

#isopentenyl diphosphate biosynthetic process, mevalonate pathway - Hmgcs1, Mvd
#detoxification of copper ion - Mt1, Mt2 (these genes are strongly linked to starvation, but direction doesnt make sense)
#geranylgeranyl diphosphate biosynthetic process - Hmgcs1, Mvd
#farnesyl diphosphate biosynthetic process, mevalonate pathway - Hmgcs1, Mvd
#zymosterol biosynthetic process - Hmgcs1, Mvd, Sqle
#cholesterol biosynthetic process via desmosterol - Hmgcs1, Mvd, Sqle
#cholesterol biosynthetic process via lathosterol - Hmgcs1, Mvd, Sqle
#retinoic acid metabolic process - Cyp2c55, Rdh16, Cyp3a11 (these are all down regulated in NY-CF the most)
#this pathway has direct links to starbation abe being downregulated during this process. RA signaling is lower in fasted mice Klyuyeva et al.2021

#Suox and Cth are related to Cysteine, which has parallels to metallothioneins (both have sulfur)


#expression ggplots for the three gene groups implicated - 1) Mt1, Mt2, 2) Hmgcs1, Mvd, Sqle, and 3) Cyp2c55, Rdh16, Cyp3a11
plot_data_mt <- map_df(c("ENSMUSG00000031765","ENSMUSG00000031762"), function(gene_id) {
  plotCounts(dds_liver,gene = gene_id,intgroup = "pop",returnData = TRUE) %>% mutate(gene = gene_id)})
plot_data_mt$gene <- factor(plot_data_mt$gene, levels = c("ENSMUSG00000031765","ENSMUSG00000031762"))

custom_labels_mt <- c("ENSMUSG00000031765" = "Mt1","ENSMUSG00000031762" = "Mt2")
ggplot(plot_data_mt, aes(x = pop, y = count, color = pop, fill = pop)) + scale_color_manual(values=paletteGroup)+
  geom_jitter(width = 0.1, size = 1) +geom_boxplot(alpha = 0.2, outlier.shape = NA) + scale_fill_manual(values=paletteGroup)+
  facet_wrap(~gene, scales = "free_y", labeller = labeller(gene = custom_labels_mt)) + labs(y = "Normalized counts", x = NULL) + 
  scale_x_discrete(labels = c("BZBZ" = "MANA","BZBZ-C" = "MANA-CF","NYNY" = "SARA", "NYNY-C" = "SARA-CF")) +
  theme_bw() + theme(legend.position = "none",
    axis.title.y = element_text(face = "bold",size = 9),
    axis.text.y = element_text(size = 8,color = "black"), 
    axis.text.x = element_text(size = 6,color = "black",face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")) 


plot_data_mvd <- map_df(c("ENSMUSG00000093930","ENSMUSG00000006517","ENSMUSG00000022351"), function(gene_id) {
  plotCounts(dds_liver,gene = gene_id,intgroup = "pop",returnData = TRUE) %>% mutate(gene = gene_id)})
plot_data_mvd$gene <- factor(plot_data_mvd$gene, levels = c("ENSMUSG00000093930","ENSMUSG00000006517","ENSMUSG00000022351"))

custom_labels_mvd <- c("ENSMUSG00000093930" = "Hmgsc1","ENSMUSG00000006517" = "Mvd", "ENSMUSG00000022351" = "Sqle")
ggplot(plot_data_mvd, aes(x = pop, y = count, color = pop, fill = pop)) + scale_color_manual(values=paletteGroup)+
  geom_jitter(width = 0.1, size = 1) +geom_boxplot(alpha = 0.2, outlier.shape = NA) + scale_fill_manual(values=paletteGroup)+
  facet_wrap(~gene, scales = "free_y", labeller = labeller(gene = custom_labels_mvd)) + labs(y = "Normalized counts", x = NULL) + 
  scale_x_discrete(labels = c("BZBZ" = "MANA","BZBZ-C" = "MANA-CF","NYNY" = "SARA", "NYNY-C" = "SARA-CF")) +
  theme_bw() + theme(legend.position = "none",
                     axis.title.y = element_text(face = "bold",size = 9),
                     axis.text.y = element_text(size = 8,color = "black"), 
                     axis.text.x = element_text(size = 6,color = "black",face = "bold"),
                     strip.text = element_text(face = "bold"),
                     strip.background = element_rect(fill = "white", color = "black")) 


plot_data_cyp <- map_df(c("ENSMUSG00000025002","ENSMUSG00000056035","ENSMUSG00000069456"), function(gene_id) {
  plotCounts(dds_liver,gene = gene_id,intgroup = "pop",returnData = TRUE) %>% mutate(gene = gene_id)})
plot_data_cyp$gene <- factor(plot_data_cyp$gene, levels = c("ENSMUSG00000025002","ENSMUSG00000056035","ENSMUSG00000069456"))

custom_labels_cyp <- c("ENSMUSG00000025002" = "Cyp2c55","ENSMUSG00000056035" = "Cyp3a11", "ENSMUSG00000069456" = "Rdh16")
ggplot(plot_data_cyp, aes(x = pop, y = count, color = pop, fill = pop)) + scale_color_manual(values=paletteGroup)+
  geom_jitter(width = 0.1, size = 1) +geom_boxplot(alpha = 0.2, outlier.shape = NA) + scale_fill_manual(values=paletteGroup)+
  facet_wrap(~gene, scales = "free_y", labeller = labeller(gene = custom_labels_cyp)) + labs(y = "Normalized counts", x = NULL) + 
  scale_x_discrete(labels = c("BZBZ" = "MANA","BZBZ-C" = "MANA-CF","NYNY" = "SARA", "NYNY-C" = "SARA-CF")) +
  theme_bw() + theme(legend.position = "none",
                     axis.title.y = element_text(face = "bold",size = 9),
                     axis.text.y = element_text(size = 8,color = "black"), 
                     axis.text.x = element_text(size = 6,color = "black",face = "bold"),
                     strip.text = element_text(face = "bold"),
                     strip.background = element_rect(fill = "white", color = "black")) 


#plot SREBf1 and SREBf2 which are critical regulators of fasting induced cholesterol changes
plot_data_srebf <- map_df(c("ENSMUSG00000020538","ENSMUSG00000022463"), function(gene_id) {
  plotCounts(dds_liver,gene = gene_id,intgroup = "pop",returnData = TRUE) %>% mutate(gene = gene_id)})
plot_data_srebf$gene <- factor(plot_data_srebf$gene, levels = c("ENSMUSG00000020538","ENSMUSG00000022463"))

custom_labels_srebf <- c("ENSMUSG00000020538" = "SREBP-1","ENSMUSG00000022463" = "SREBP-2")
ggplot(plot_data_srebf, aes(x = pop, y = count, color = pop, fill = pop)) + scale_color_manual(values=paletteGroup)+
  geom_jitter(width = 0.1, size = 1) +geom_boxplot(alpha = 0.2, outlier.shape = NA) + scale_fill_manual(values=paletteGroup)+
  facet_wrap(~gene, scales = "free_y", labeller = labeller(gene = custom_labels_srebf)) + labs(y = "Normalized counts", x = NULL) + 
  scale_x_discrete(labels = c("BZBZ" = "MANA","BZBZ-C" = "MANA-CF","NYNY" = "SARA", "NYNY-C" = "SARA-CF")) +
  theme_bw() + theme(legend.position = "none",
                     axis.title.y = element_text(face = "bold",size = 9),
                     axis.text.y = element_text(size = 7,color = "black"), 
                     axis.text.x = element_text(size = 8,color = "black",face = "bold"),
                     strip.text = element_text(face = "bold"),
                     strip.background = element_rect(fill = "white", color = "black")) 



#plot the GO
#biosynthetic process is abbreviated BP
hub_GOterms <- read.csv("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/liverMatHubGO.csv")

ggplot(data=hub_GOterms, aes(x=reorder(GO.biological.process.complete, -log10(FDR)), y=-log10(FDR))) + 
  scale_x_discrete(position = "top") + scale_y_reverse(limits=c(1.3,3.2)) +geom_segment(aes(yend=1.3, xend = GO.biological.process.complete)) +
  geom_point(aes(size=Fold.Enrichment)) + coord_flip() + labs(size="Fold \nEnrichment") +
  theme_classic() +  ylab("-Log10 FDR") +theme(legend.position = 'bottom', 
                                                    axis.title.x = element_text(face="bold", size =12),
                                                    axis.text.x = element_text(size =10, colour = 'black'), 
                                                    axis.title.y = element_blank(),
                                                    axis.text.y = element_text(size =9, colour = 'black'))






#### HEAT MAP ####
library(pheatmap)
library(RColorBrewer)
library(viridis)

#renormalize
normalized_dds_liver <- assay(dds_liver_norm) 

#subset dds for all module and hub gene heatmaps
mat_blueMod_liver <- normalized_dds_liver[rownames(normalized_dds_liver) %in% blueModule_genes, ]
mat_hub_liver <- normalized_dds_liver[rownames(normalized_dds_liver) %in% hubCandidates_matLiver, ]

#write function to replace gene names with common symbols
replace_ids_with_symbols <- function(mat_sub, gene_map) {
  map_order <- match(rownames(mat_sub), gene_map$gene)
  rownames(mat_sub) <- gene_map$gene_name[map_order]
  
  # remove duplicates
  mat_sub <- mat_sub[!duplicated(rownames(mat_sub)), ]
  
  return(mat_sub)
}

#replace gene names with common symbols
mat_blueMod_liver <- replace_ids_with_symbols(mat_blueMod_liver, modGenes_named)
mat_hub_liver <- replace_ids_with_symbols(mat_hub_liver, hubGenes_named)

#scale for visualization
mat_blueMod_liver_scaled <- t(scale(t(mat_blueMod_liver)))
mat_hub_liver_scaled <- t(scale(t(mat_hub_liver)))

#add annotation column for groups
annotation_col_blueMod <- data.frame(Group = colData(dds_liver)$pop)
rownames(annotation_col_blueMod) <- colnames(mat_blueMod_liver_scaled)
annotation_col_hub <- data.frame(Group = colData(dds_liver)$pop)
rownames(annotation_col_hub) <- colnames(mat_hub_liver_scaled)


#plot - just a basic version
ann_colors = list(Group = c("BZBZ" = "#ce8c14","BZBZ-C" = "#f2cd8b","NYNY" = "#006d92","NYNY-C" = "#7fbcd1"))

pheatmap(mat_blueMod_liver_scaled, annotation_col = annotation_col_blueMod, cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, color = inferno(100), border_color = NA,
         annotation_colors = ann_colors)

pheatmap(mat_hub_liver_scaled, annotation_col = annotation_col_hub, cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, color = inferno(100), border_color = NA,
         annotation_colors = ann_colors)








