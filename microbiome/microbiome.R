##Fecal microbiome analyses from CF and non-CF SARA and MANA mice

#we have species level and genus level data - but will focus on genus level as it is more reliable
#the species level is a mess...didnt really work


#set up
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ggpubr)
library(RColorBrewer)
library(paletteer)
library(ggrepel)

setwd("~/Documents/Berkeley/Research/Maternal Effects/Microbiomes")
genus_counts <- read.table("genus_abundance_matrix.tsv", sep="\t", header=TRUE) #211 genera represented
rownames(genus_counts) <- genus_counts[,1]
genus_counts <- genus_counts[,-1]

#filter genera for a total count across samples greater than 10 and being prevalent in at least 4 samples (same size as an experimental group)
genus_counts_filt <- genus_counts %>% filter(rowSums(genus_counts >= 10) >= 4) #37 genera represented

#construct metadata table
metadata <- read.csv("microbiome_metadata.csv")

#filter data and metadata for pups only
pup_metadata <- subset(metadata, age == "pup")
pup_genus_counts <- genus_counts_filt[, pup_metadata$id, drop=FALSE]




#### PCAs / Distributions ####

#transpose counts
genus_counts_filt_t <- t(genus_counts_filt)
pup_genus_counts_t <- t(pup_genus_counts)


#normalize counts to be relative abundance - then log transform
all_rel <- genus_counts_filt_t / rowSums(genus_counts_filt_t)
pup_rel <- pup_genus_counts_t / rowSums(pup_genus_counts_t)

all_log <- log10(all_rel + 1e-6)
pup_log <- log10(pup_rel + 1e-6)


#run the PCA and structure the dataframe
all_pca_res <- prcomp(all_log, center = TRUE, scale. = TRUE)
pup_pca_res <- prcomp(pup_log, center = TRUE, scale. = TRUE)

all_pca_df <- as.data.frame(all_pca_res$x)
all_pca_df$id <- rownames(all_pca_df)
all_pca_df <- merge(all_pca_df, metadata, by = "id")

pup_pca_df <- as.data.frame(pup_pca_res$x)
pup_pca_df$id <- rownames(pup_pca_df)
pup_pca_df <- merge(pup_pca_df, pup_metadata, by = "id")


#calculation percent variance explained
all_percent_var <- round(100 * (all_pca_res$sdev^2 / sum(all_pca_res$sdev^2)), 1)
pup_percent_var <- round(100 * (pup_pca_res$sdev^2 / sum(pup_pca_res$sdev^2)), 1)


#plots!
paletteAll <- c("#ce8c14","#f2cd8b","grey60","gray80","black","gray40", "#006d92", "#7fbcd1")
paletteGroup <- c("#ce8c14","#f2cd8b", "#006d92", "#7fbcd1")

#all samples
ggplot(all_pca_df, aes(PC1, PC2, color = group, shape = age)) +geom_point(size = 4) +scale_color_manual(values=paletteAll)+
  theme_bw() + xlab(paste0("PC1 (", all_percent_var[1], "%)")) + ylab(paste0("PC2 (", all_percent_var[2], "%)"))

#pups only - definitely cleaner grouping than the Epi2me based analysis - add in driving taxa (sig diff abundance between rearing OR pop - done below)
#add the significant taxa
loadings <- as.data.frame(pup_pca_res$rotation)
loadings$Taxon <- rownames(loadings)
sig_taxa <- c("Lactobacillus","Klebsiella","HT002","Anaerostipes","Tyzzerella","Anaeroplasma", "Blautia", "Faecalimonas","Butyribacter","Bacteroides")
loadings_sig <- loadings[loadings$Taxon %in% sig_taxa, ]

#scale the loadings
loadings_sig$PC1 <- loadings_sig$PC1 * 10
loadings_sig$PC2 <- loadings_sig$PC2 * 10

#plot
ggplot(pup_pca_df, aes(PC1, PC2, color = group, fill=group)) +scale_color_manual(values = paletteGroup) +scale_fill_manual(values = paletteGroup) +
  xlab(paste0("PC1 (", pup_percent_var[1], "%)")) +ylab(paste0("PC2 (", pup_percent_var[2], "%)")) +
  geom_segment(data = loadings_sig,aes(x = 0, y = 0, xend = PC1, yend = PC2),inherit.aes = FALSE,arrow = arrow(length = unit(0.2, "cm")),color = "grey70") +# arrows
  geom_text_repel(data = loadings_sig,aes(x = PC1, y = PC2, label = Taxon),inherit.aes = FALSE,size = 2.5,color = "grey30") +# labels
  geom_point(size = 4, alpha=0.85) + theme_bw() +  theme(legend.position = 'none', 
                       axis.title.x = element_text(face="bold", size =10),
                       axis.text.x = element_text(size =8, colour = 'black'), 
                       axis.title.y = element_text(face="bold", size =10),
                       axis.text.y = element_text(size =8, colour = 'black'))


#PERMANOVA - statistics on clustering / what variables are explaining variation
library(vegan)
library(ape)

# Bray-Curtis distance (better for microbiome than PCA distance)
dist_mat <- vegdist(pup_rel, method = "bray")

adonis2(dist_mat ~ pop + rearingDam, data = pup_metadata, by="margin")
#p value for pop is 0.022 and explains 13.6% of var, p value for rearingDam is 0.004 and explains 45% of var
adonis2(dist_mat ~ group, data = pup_metadata, by="margin")
#p value for group is 0.001 and explains 65% of var

#pairwise tests from above - PERMANOVAs
library(pairwiseAdonis)
pairwise.adonis2( dist_mat ~ group ,data = pup_metadata,p.adjust.m = "BH")
#NS:only BZvNY-CF

#check dispersion and if that is driving our differences
disp <- betadisper(dist_mat, pup_metadata$group)
anova(disp) #p>0.05
permutest(disp) #p>0.05



#make a PCoA plot of bray-curtis dissimilarity
pcoa_pups <- pcoa(dist_mat)
pcoa_pups_df <- as.data.frame(pcoa_pups$vectors[,1:2])
pcoa_pups_df$id <- rownames(pcoa_pups_df)
pcoa_pups_df <- merge(pcoa_pups_df, pup_metadata, by = "id")

#plot
ggplot(pcoa_pups_df, aes(Axis.1, Axis.2, color = group, shape = fosterStatus))+geom_point(size = 4) + scale_color_manual(values=paletteGroup)+theme_bw()

#boxplot of bray-curtis dissimilarilty
bc_mat <- as.matrix(dist_mat)
bc_long <- as.data.frame(as.table(bc_mat))
colnames(bc_long) <- c("sample1", "sample2", "BrayCurtis")
bc_long <- bc_long %>%filter(as.character(sample1) < as.character(sample2))

#add metadata
meta_small <- pup_metadata %>% select(id, group, pop, rearingDam, fosterStatus)
bc_long <- bc_long %>%filter(as.character(sample1) < as.character(sample2)) %>%
  left_join(meta_small, by = c("sample1" = "id")) %>%left_join(meta_small, by = c("sample2" = "id"), suffix = c("1", "2"))

#calc distances to plot
bc_between <- bc_long %>%filter(group1 != group2) %>%mutate(pair = paste(group1, group2, sep = " vs "))

#plot! - palette is defined below
ggplot(bc_between, aes(x = reorder(pair, -BrayCurtis), y = BrayCurtis, fill=reorder(pair, -BrayCurtis))) +geom_boxplot(alpha=0.7, width=0.65) +
  geom_jitter(aes(color=reorder(pair, -BrayCurtis)),width = 0.1, alpha = 0.6) + ylim(0.16, 1)+
  theme_bw() +ylab("Bray-Curtis dissimilarity") +xlab("") + scale_color_manual(values = new_pal) + scale_fill_manual(values = new_pal)+
  theme(legend.position = 'none',
              axis.text.x = element_text(face="bold", size =8, colour = 'black',angle = 45, hjust = 1), 
              axis.title.y = element_text(face="bold", size =10), 
              axis.title.x = element_blank(), 
              axis.text.y = element_text(size =8, colour = 'black'))





#### Differential abundance ####

#use DESeq to quantify differential abundace - both between rearing status + pop  and separate group comparisons
dds_pupMicro_rearing <- DESeqDataSetFromMatrix(countData = pup_genus_counts, colData = pup_metadata, design = ~ pop + rearingDam)
dds_pupMicro_group <- DESeqDataSetFromMatrix(countData = pup_genus_counts, colData = pup_metadata, design = ~ group)

#run deseq
dds_pupMicro_rearing <- DESeq(dds_pupMicro_rearing)
dds_pupMicro_group <- DESeq(dds_pupMicro_group)

#get different results dataframes
#overall rearing and pop
res_rearing <- as.data.frame(results(dds_pupMicro_rearing, contrast = c("rearingDam", "NY", "BZ"))) # neg LFC = higher in BZ reared
res_pop <- as.data.frame(results(dds_pupMicro_rearing, contrast = c("pop", "NY", "BZ"))) # neg LFC = higher in BZ the pop

#specific groups
res_NYvsBZ <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "NYNY", "BZBZ"))) # neg LFC = higher in BZ 
res_NYCFvsBZCF <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "NYNY-CF", "BZBZ-CF"))) # neg LFC = higher in BZ-CF
res_NYVNYCF <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "NYNY", "NYNY-CF"))) # neg LFC = higher in NY-CF
res_NYVBZCF <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "NYNY", "BZBZ-CF"))) # neg LFC = higher in BZ-CF
res_BZCFVBZ <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "BZBZ-CF", "BZBZ"))) # neg LFC = higher in BZ 
res_NYCFvsBZ <- as.data.frame(results(dds_pupMicro_group, contrast = c("group", "NYNY-CF", "BZBZ"))) # neg LFC = higher in BZ 

#Lactobacillus is the highest differential in rearing and NYvNYCF - it just isnt present in NY-CF or BZ pups - similarly it is not in BZ dams but is high in NY dams
ggplot(plotCounts(dds_pupMicro_group, gene = "Lactobacillus", intgroup = "group", returnData=TRUE), aes(x=group, y=count, color=group)) +scale_color_manual(values=paletteGroup)+
  geom_point(position=position_jitter(w=0.1,h=0), size=2) + ggtitle("Lactobacillus abundance") + theme_bw() #related to cholesterol




#### overall abundance stacked plots ####

#format dataframes and integrate metadata
all_rel_long <- as.data.frame(all_rel) %>% rownames_to_column("id") %>% pivot_longer(-id, names_to = "Taxon", values_to = "Abundance")
pup_rel_long <- as.data.frame(pup_rel) %>% rownames_to_column("id") %>% pivot_longer(-id, names_to = "Taxon", values_to = "Abundance")

plot_all_rel <- as.data.frame(left_join(all_rel_long, metadata, by = "id"))
plot_pup_rel <- as.data.frame(left_join(pup_rel_long, pup_metadata, by = "id"))

#focus in on just the top 10 most abundant taxa
all_top_taxa <- plot_all_rel %>%group_by(Taxon) %>%summarise(total = sum(Abundance), .groups = "drop") %>%arrange(desc(total)) %>%head(14) %>%pull(Taxon)
pup_top_taxa <- plot_pup_rel %>%group_by(Taxon) %>%summarise(total = sum(Abundance), .groups = "drop") %>%arrange(desc(total)) %>%head(14) %>%pull(Taxon)

#collapse the remainder into "other"
plot_all_rel2 <- plot_all_rel %>% mutate(Taxon2 = ifelse(Taxon %in% all_top_taxa, Taxon, "Other")) %>%
  group_by(id, group, plot_group,pop, rearingDam, fosterStatus, sex, Taxon2) %>% summarise(Abundance = sum(Abundance), .groups = "drop")

plot_pup_rel2 <- plot_pup_rel %>% mutate(Taxon2 = ifelse(Taxon %in% pup_top_taxa, Taxon, "Other")) %>%
  group_by(id, group, pop, rearingDam, fosterStatus, sex, Taxon2) %>% summarise(Abundance = sum(Abundance), .groups = "drop")

#order samples
plot_all_rel2$id <- factor(plot_all_rel2$id,levels = unique(plot_all_rel2$id[order(plot_all_rel2$group)]))
plot_pup_rel2$id <- factor(plot_pup_rel2$id,levels = unique(plot_pup_rel2$id[order(plot_pup_rel2$group)]))

#order the legend
plot_all_rel2$Taxon2 <- factor(plot_all_rel2$Taxon2,levels = c(all_top_taxa, "Other"))
plot_pup_rel2$Taxon2 <- factor(plot_pup_rel2$Taxon2,levels = c(pup_top_taxa, "Other"))

#reorder gourps for all plots
plot_all_rel2$plot_group <- factor(plot_all_rel2$plot_group,levels = c("BZdam", "BZBZ", "BZBZ-CF","NYdam", "NYNY", "NYNY-CF"))

#plots!
my_pal <- as.vector(paletteer_d("MetBrewer::Signac", direction=-1))
new_pal <- c(my_pal, "grey")

ggplot(plot_all_rel2, aes(x = id, y = Abundance, fill = Taxon2)) +geom_bar(stat = "identity", width = 1,alpha=0.8) +
  facet_grid(~ plot_group, scales = "free_x", space = "free_x") +scale_fill_manual(values = new_pal)+
  scale_y_continuous(expand = c(0.0005, 0.01))+theme_classic()+ theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank(),strip.text = element_text(size = 11, face = "bold"),
                         axis.ticks = element_blank(),strip.background = element_blank(),
                         panel.spacing.x = unit(0.01, "lines"), legend.title = element_blank())

ggplot(plot_pup_rel2, aes(x = id, y = Abundance, fill = Taxon2)) +geom_bar(stat = "identity", width = 1,alpha=0.8) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +scale_fill_manual(values = new_pal)+
  scale_y_continuous(expand = c(0.0005, 0.01))+theme_classic() +theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank(),strip.text = element_text(size = 11, face = "bold"),
                         axis.ticks = element_blank(),strip.background = element_blank(),
                         panel.spacing.x = unit(0.01, "lines"), legend.title = element_blank())




#### Diversity analyses ####
# Shannon
metadata$Shannon <- diversity(all_rel, index = "shannon")
# Simpson
metadata$Simpson <- diversity(all_rel, index = "simpson")
# Richness (presence/absence)
metadata$Richness <- specnumber(all_rel > 0)

#subset pups here
pups_div <- subset(metadata, age == "pup")

#statistical tests - by group and rearing dam
anova(lm(Shannon ~ pop + rearingDam, data = pups_div)) #NS
anova(lm(Shannon ~ group, data = pups_div)) #NS
#Tukey post-hoc test
fit_pupShannon <- aov(Shannon ~ group, data = pups_div)
TukeyHSD(fit_pupShannon) 
#plot
ggplot(pups_div, aes(x = group, y = Shannon, color=group)) +geom_boxplot() +geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c("BZBZ", "BZBZ-CF"), c("NYNY", "NYNY-CF"), c("BZBZ", "NYNY")), label = "p.format", size=3)+ 
  scale_color_manual(values=paletteGroup)+theme_bw() #pvalues are a little higher in the correct statistical test above

anova(lm(Simpson ~ pop + rearingDam, data = pups_div)) #NS
anova(lm(Simpson ~ group, data = pups_div)) #NS
#Tukey post-hoc test
fit_pupSimpson <- aov(Simpson ~ group, data = pups_div)
TukeyHSD(fit_pupSimpson) 
#plot
ggplot(pups_div, aes(x = group, y = Simpson, color=group)) +geom_boxplot() +geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c("BZBZ", "BZBZ-CF"), c("NYNY", "NYNY-CF"), c("BZBZ", "NYNY")), label = "p.format", size=3)+ 
  scale_color_manual(values=paletteGroup)+theme_bw() #pvalues are a little higher in the correct statistical test above

anova(lm(Richness ~ pop + rearingDam, data = pups_div)) #p=0.01457
anova(lm(Richness ~ group, data = pups_div)) #p=0.01117
#Tukey post-hoc test
fit_pupRich <- aov(Richness ~ group, data = pups_div)
TukeyHSD(fit_pupRich) 
#the only different one is NYNY-BZBZ (NY more rich)...not any CF, those are 0.06 p values though (NYNY-CF-NYNY and NYNY-BZBZ-CF)

#plot this
ggplot(pups_div, aes(x = group, y = Richness, color=group)) +geom_boxplot() +geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c("BZBZ", "BZBZ-CF"), c("NYNY", "NYNY-CF"), c("BZBZ", "NYNY")), label = "p.format", size=3, method = "t.test")+ 
  scale_color_manual(values=paletteGroup)+theme_bw() #pvalues are a little higher in the correct statistical test above


#test in dams v pups
anova(lm(Shannon ~ age, data = metadata)) #NS
anova(lm(Shannon ~ plot_group, data = metadata)) #NS

anova(lm(Simpson ~ age, data = metadata)) #p=0.03675
anova(lm(Simpson ~ plot_group, data = metadata)) #NS

anova(lm(Richness ~ age, data = metadata)) #p=0.0003386
anova(lm(Richness ~ plot_group, data = metadata)) #p=2.038e-05
#Tukey post-hoc test
fit_allRich <- aov(Richness ~ plot_group, data = metadata)
TukeyHSD(fit_allRich) 
#dam to pup are mostly all different and NYdam is diff from BZdam (NY more rich)





#### F:B ratio ####
#the Firmicutes/Bacteroidetes (F/B) ratio is commonly cited as a marker of obesity/weight gain
#these are different phyla so we will have to group our represented genera to assess this

firmicutes_genera <- c("Marvinbryantia", "Staphylococcus", "Enterococcus","Acutalibacter", "Paludicola", "Lachnoclostridium",
  "Tyzzerella", "Blautia", "Agathobaculum", "Anaerotruncus", "Limosilactobacillus", "Oscillibacter","Lactobacillus", 
  "Acetitomaculum", "Faecalimonas", "Roseburia", "Colidextribacter", "Ligilactobacillus","Streptococcus", "Anaerostipes", "Butyribacter")
bacteroidetes_genera <- c("Muribaculum", "Bacteroides", "Rikenella","Alistipes", "Odoribacter")

#count and calculate ratio
firm_counts <- colSums(genus_counts_filt[firmicutes_genera, , drop = FALSE])
bact_counts <- colSums(genus_counts_filt[bacteroidetes_genera, , drop = FALSE])

fb_ratio <- firm_counts / bact_counts

#integrate into metadata
metadata$FB_ratio <- fb_ratio[match(metadata$id, names(fb_ratio))]

#statistically test
anova(lm(FB_ratio ~ pop + rearingDam, data = subset(metadata, age=="pup"))) #both are sig, pop more so
anova(lm(FB_ratio ~ group, data = subset(metadata, age=="pup"))) #p=0.003365
#Tukey post-hoc test
fit_pupFB <- aov(Richness ~ group, data = subset(metadata, age=="pup"))
TukeyHSD(fit_pupFB) #higher in our NYNY mice, but this is statistically not different from NYNY-CF (BUT p=0.069), but it is significantly higher than BZBZ

#plot the F:B ratio per group
ggplot(metadata, aes(x = group, y = FB_ratio, color=group)) + geom_boxplot() +geom_jitter(width = 0.1) +theme_bw()

ggplot(subset(metadata, age == "pup"), aes(x = group, y = FB_ratio, color=group)) +geom_boxplot() +geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c("BZBZ", "BZBZ-CF"), c("NYNY", "NYNY-CF"), c("BZBZ", "NYNY"),c("BZBZ", "NYNY-CF")), label = "p.format", size=3)+ 
  scale_color_manual(values=paletteGroup)+theme_bw()
#muuuuch higher in our NYNY mice, but this is statistically not different from NYNY-CF (BUT p=0.057), but it is significantly higher than BZBZ





#similarly to Richness, F:B is highest in NYNY mice, which is significantly higher than BZBZ, and substantially but not significantly higher than the NYNY-CF counterpart (p~0.06 for both)
