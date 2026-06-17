#contains code to define regulatory categories (cis v trans) in F1 ASE data, liver only, sexes grouped

#contains code to make plots in Figures 2B, C, and E, and Supplementary Figure 9

#result dataframe names for each test type
#P: parent alone wald's test: parents_NYBZ_liver, parents_BZNY_liver
#H: hybrid alone wald's test: NYBZ_ASE_liver, BZNY_ASE_liver
#T: trans component LRT test: F1andParent_NYBZ_liver, F1andParent_BZNY_liver

#Main dataframes are merge_results_NYBZ_liver and merge_results_BZNY_liver



############################################################## #
##### Set up ####
############################################################## #

#Navigate to correct folder and load ASE data 
setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_ASE_comps.RData")

library(ggplot2)
library(ggrastr)

#mapping bias histograms
NYBZ.liverMappingBias <- NYBZ_liver_countData
NYBZ.liverMappingBias$ratio <- (rowSums(NYBZ.liverMappingBias[c(2:13)])) / (rowSums(NYBZ.liverMappingBias[c(2:25)])) #2 through 13 are the BZ reads, and 2 - 25 is all reads
NYBZ.liverMappingBias <- na.omit(NYBZ.liverMappingBias)

BZNY.liverMappingBias <- BZNY_liver_countData
BZNY.liverMappingBias$ratio <- (rowSums(BZNY.liverMappingBias[c(2:13)])) / (rowSums(BZNY.liverMappingBias[c(2:25)])) #2 through 13 are the BZ reads, and 2 - 25 is all reads
BZNY.liverMappingBias <- na.omit(BZNY.liverMappingBias)

#plots
ggplot(NYBZ.liverMappingBias, aes(x=ratio)) + 
  geom_histogram(color="black", fill="grey89", binwidth = 0.01) + theme_classic() + labs(title="SARA x MANA liver",x="REF allele / (REF + ALT allele)", y = "Count") + 
  theme(axis.line.x = element_line(colour = 'black', size=0.7, linetype='solid'), axis.title.x = element_text(face="bold", size =12),
        axis.text.x = element_text(size =9, colour = 'black'), axis.line.y = element_line(colour = 'black', size=0.7, linetype='solid'),
        axis.title.y = element_text(face="bold", size =12),axis.text.y = element_text(size =9, colour = 'black'))

ggplot(BZNY.liverMappingBias, aes(x=ratio)) + 
  geom_histogram(color="black", fill="grey89", binwidth = 0.01) + theme_classic() + labs(title="MANA x SARA liver",x="REF allele / (REF + ALT allele)", y = "Count") + 
  theme(axis.line.x = element_line(colour = 'black', size=0.7, linetype='solid'), axis.title.x = element_text(face="bold", size =12),
        axis.text.x = element_text(size =9, colour = 'black'), axis.line.y = element_line(colour = 'black', size=0.7, linetype='solid'),
        axis.title.y = element_text(face="bold", size =12),axis.text.y = element_text(size =9, colour = 'black'))


############################################################## #
##### Defining genes falling into each regulatory category ####
############################################################## #

#cis only: sig in P and H, not sig in T
#trans only: sig in P, not sig in H, sig in T
#cis + trans: sig in P, H and T. sign of log2FC in P and H is the SAME direction, same direction
#cis x trans: sig in P, H and T. sign of log2FC in P and H is the OPPOSITE direction, opposing 
#compensatory: sig in H, not sig in P, sig in T. no expression divergence due to compensation
#conserved: not sig in H, P, or T. no expression divergence
#ambiguous: the left over potential patterns.


##### CIS ONLY: sig in H and P, not sig in T

#NYBZ - 1305 of 7679 total total
cisOnly_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj > 0.05))
#BZNY - 1213 of 7471 total total
cisOnly_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj > 0.05))


##### TRANS ONLY: sig in P, not sig in H, sig in T

#NYBZ - 761 of 7679 total total
transOnly_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj < 0.05))
#BZNY - 738 of 7471 total total
transOnly_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj > 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj < 0.05))


##### CIS + TRANS: sig in P, H and T. sign of log2FC in P and H is the SAME direction, same direction 

#NYBZ - 404 of 7679
cisANDtrans_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj < 0.05) &
                             (((merge_results_NYBZ_liver$parent.NYvBZ_Log2FC < 0) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_Log2FC < 0)) | ((merge_results_NYBZ_liver$parent.NYvBZ_Log2FC > 0) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_Log2FC > 0))))
#BZNY - 432 of 7471
cisANDtrans_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj < 0.05) &
                             (((merge_results_BZNY_liver$parent.NYvBZ_Log2FC < 0) & (merge_results_BZNY_liver$BZNY.NYvBZ_Log2FC < 0)) | ((merge_results_BZNY_liver$parent.NYvBZ_Log2FC > 0) & (merge_results_BZNY_liver$BZNY.NYvBZ_Log2FC > 0))))


##### CIS x TRANS: sig in P, H and T. sign of log2FC in P and H is the OPPOSITE direction, opposing 

#NYBZ - 177 of 7679
cisXtrans_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj < 0.05) &
                           (((merge_results_NYBZ_liver$parent.NYvBZ_Log2FC < 0) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_Log2FC > 0)) | ((merge_results_NYBZ_liver$parent.NYvBZ_Log2FC > 0) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_Log2FC < 0))))
#BZNY - 187 of 7471
cisXtrans_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj < 0.05) &
                           (((merge_results_BZNY_liver$parent.NYvBZ_Log2FC < 0) & (merge_results_BZNY_liver$BZNY.NYvBZ_Log2FC > 0)) | ((merge_results_BZNY_liver$parent.NYvBZ_Log2FC > 0) & (merge_results_BZNY_liver$BZNY.NYvBZ_Log2FC < 0))))


####### COMPENSATORY: sig in H, not sig in P, sig in T. no expression divergence due to compensation 

#NYBZ - 322 of 7679
compensatory_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj < 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj < 0.05))
#BZNY - 337 of 7471
compensatory_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj > 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj < 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj < 0.05))


##### CONSERVED: not sig in H, P, or T. no expression divergence 

#NYBZ - 2308 of 7679
conserved_NYBZ_liver <- subset(merge_results_NYBZ_liver, (merge_results_NYBZ_liver$parent.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver$NYBZ.NYvBZ_padj > 0.05) & (merge_results_NYBZ_liver$int.NYvBZ_padj > 0.05))
#BZNY - 2269 of 7471
conserved_BZNY_liver <- subset(merge_results_BZNY_liver, (merge_results_BZNY_liver$parent.NYvBZ_padj > 0.05) & (merge_results_BZNY_liver$BZNY.NYvBZ_padj > 0.05) & (merge_results_BZNY_liver$int.NYvBZ_padj > 0.05))




#### make combo dataset 
#has genes with their mode of regulation specified 

#first change rownames to gene names and add a column in each dataset for mode of regulation
#cis only
cisOnly_NYBZ_liver$regulation_mode <- "cis only"
cisOnly_BZNY_liver$regulation_mode <- "cis only"
#trans only
transOnly_NYBZ_liver$regulation_mode <- "trans only"
transOnly_BZNY_liver$regulation_mode <- "trans only"
#cis+trans
cisANDtrans_NYBZ_liver$regulation_mode <- "cis + trans"
cisANDtrans_BZNY_liver$regulation_mode <- "cis + trans"
#cisXtrans
cisXtrans_NYBZ_liver$regulation_mode <- "cis x trans"
cisXtrans_BZNY_liver$regulation_mode <- "cis x trans"
#compensatory
compensatory_NYBZ_liver$regulation_mode <- "compensatory"
compensatory_BZNY_liver$regulation_mode <- "compensatory"
#conserved
conserved_NYBZ_liver$regulation_mode <- "conserved"
conserved_BZNY_liver$regulation_mode <- "conserved"


#bind the above datasets to have one for each comparison (NYvBZ.liver, NYvBZ.BAT. CANvBZ.liver, CANvBZ.BAT)
NYBZ.liver_RegulationModes <- rbind(cisOnly_NYBZ_liver, transOnly_NYBZ_liver, cisANDtrans_NYBZ_liver, cisXtrans_NYBZ_liver, compensatory_NYBZ_liver, conserved_NYBZ_liver) #5277 of 7679
BZNY.liver_RegulationModes <- rbind(cisOnly_BZNY_liver, transOnly_BZNY_liver, cisANDtrans_BZNY_liver, cisXtrans_BZNY_liver, compensatory_BZNY_liver, conserved_BZNY_liver) #5176 0f 7471


#write these datasets to my computer for later access 
write.csv(NYBZ.liver_RegulationModes, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/NYBZ.liver_RegulationModes.csv", row.names = TRUE)
write.csv(BZNY.liver_RegulationModes, "~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/BZNY.liver_RegulationModes.csv", row.names = TRUE)




#### Plotting ####

#scatterplots of regulatory types for NYBZ and BZNY separately
#NYBZ
ggplot(merge_results_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC)) + 
  rasterize(geom_point(size = 0.6, color="#bab0ac", alpha=0.7), dpi = 150) + xlim(-3.5,3.5) + ylim(-3.5,3.5) +
  rasterize(geom_point(data=compensatory_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC), fill="#3399ff", color="#236bb2", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisANDtrans_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC), fill="#ffae42", color="#cc8b34", size=1.5, shape=21,  alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisXtrans_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC), fill="#03c04a",color= "#028633", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=transOnly_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC), fill="#9966cd",color= "#7a51a4", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisOnly_NYBZ_liver, aes(x=parent.NYvBZ_Log2FC, y=NYBZ.NYvBZ_Log2FC), fill="#f95738",color="#c7452c", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  theme_bw() + ylab(bquote(bold(SARA/MANA~allele~log[2]~Fold~Change~(F1s)))) + xlab(bquote(bold(SARA/MANA~log[2]~Fold~Change~(parents)))) + theme(
    axis.text.x = element_text( size =9, colour = 'black'), 
    axis.title.y = element_text(face="bold", size =11),
    axis.title.x = element_text(face="bold", size =11),
    axis.text.y = element_text(size =9, colour = 'black'))

#BZNY
ggplot(merge_results_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC)) + 
  rasterize(geom_point(size = 0.6, color="#bab0ac", alpha=0.7), dpi = 150) + xlim(-3.5,3.5) + ylim(-3.5,3.5) +
  rasterize(geom_point(data=compensatory_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC), fill="#3399ff", color="#236bb2", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisANDtrans_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC), fill="#ffae42", color="#cc8b34", size=1.5, shape=21,  alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisXtrans_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC), fill="#03c04a",color= "#028633", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=transOnly_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC), fill="#9966cd",color= "#7a51a4", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  rasterize(geom_point(data=cisOnly_BZNY_liver, aes(x=parent.NYvBZ_Log2FC, y=BZNY.NYvBZ_Log2FC), fill="#f95738",color="#c7452c", size=1.5, shape=21, alpha=0.5), dpi = 150) +
  theme_bw() + ylab(bquote(bold(SARA/MANA~allele~log[2]~Fold~Change~(F1s)))) + xlab(bquote(bold(SARA/MANA~log[2]~Fold~Change~(parents)))) + theme(
    axis.text.x = element_text( size =9, colour = 'black'), 
    axis.title.y = element_text(face="bold", size =11),
    axis.title.x = element_text(face="bold", size =11),
    axis.text.y = element_text(size =9, colour = 'black'))





#### Regulation of maternal and genetic sets ####

#mat_either_genes_liver and gen_genes_liver object is made in LIVERmat_vs_gen.R

setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/Transferred Files/")
#load regulatory modes spreadsheet - contains the cis only regulated genes
NYBZ.liver_RegulationModes <- read.csv("NYBZ.liver_RegulationModes.csv")
BZNY.liver_RegulationModes <- read.csv("BZNY.liver_RegulationModes.csv")

#create dataset that has reg modes of mat genes
NYBZregModes_mat_either_genes_liver <- na.omit(merge(NYBZ.liver_RegulationModes,mat_either_genes_liver,by="gene", all.y = TRUE)) #74, was 55 have reg modes
BZNYregModes_mat_either_genes_liver <- na.omit(merge(BZNY.liver_RegulationModes,mat_either_genes_liver,by="gene", all.y = TRUE)) #74, was 56 have reg modes

#create dataset that has reg modes of genetic genes
NYBZregModes_gen_genes_liver <- na.omit(merge(NYBZ.liver_RegulationModes,gen_genes_liver,by="gene", all.y = TRUE)) #677 have reg modes
BZNYregModes_gen_genes_liver <- na.omit(merge(BZNY.liver_RegulationModes,gen_genes_liver,by="gene", all.y = TRUE)) #647 have reg modes

#Break down of regulatory types in maternal sets
NYBZregModes_mat_either_liver_regProps <- as.data.frame(prop.table(table(NYBZregModes_mat_either_genes_liver$regulation_mode)) * 100)
#cis = (39.2%), trans = (31.1%), cis + trans = (21.6%), cis x trans = (8.1%)
BZNYregModes_mat_either_liver_regProps <- as.data.frame(prop.table(table(BZNYregModes_mat_either_genes_liver$regulation_mode)) * 100)
#cis = (32.4%), trans = (28.4%), cis + trans = (29.7%), cis x trans = (9.5%)

#Break down of regulatory types in GENETIC sets
NYBZregModes_gen_liver_regProps <- as.data.frame(prop.table(table(NYBZregModes_gen_genes_liver$regulation_mode)) * 100)
#cis = (53.3%), trans =  (20.3%), cis + trans =  24.2%), cis x trans =  (2.2%)
BZNYregModes_gen_liver_regProps <- as.data.frame(prop.table(table(BZNYregModes_gen_genes_liver$regulation_mode)) * 100)
#cis = (52.8%), trans = (17.9%), cis + trans = (27.2%), cis x trans = (2.1%)

#overall NYBZ breakdown (out of the four categories of divergent parent genes)
#cis = 1290/2632 (49.0%), trans = 761/2632 (28.9%), cis + trans = 406/2632 (15.4%), cis x trans = 175/2632 (6.6%)


#plot it!
library(ggalluvial)
library(scales)

#separate alluvial plots for NYBZ and BZNY of liver mat and genetic regulatory breakdown

regPercents_NYBZ <- data.frame(comparison = c("Maternal","Maternal", "Maternal", "Maternal", "Genetic", "Genetic","Genetic","Genetic"),
                                reg_type = c("cis only", "trans only", "cis + trans", "cis x trans","cis only", "trans only", "cis + trans", "cis x trans"),
                                reg_percent = c(39.2,31.1,21.6,8.1,   53.3,20.3,24.4,2.2))

ggplot(regPercents_NYBZ, aes(fill=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")), y=reg_percent, x=comparison, 
                                    alluvium=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")), 
                                    stratum=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")))) + 
  geom_flow(alpha=0.6,width = 0.6)+  geom_stratum(color = "transparent",alpha=0.95, width = 0.6)+
  geom_text(stat = "stratum",aes(label = percent(after_stat(prop), accuracy = .1)), size = 3) + scale_fill_manual("Regulation",values = c("#f95738", "#9966cd","#ffae42", "#03c04a")) +
  scale_x_discrete(limits=c("Genetic","Maternal")) +theme_classic() + ylab("Proportion of Regulatory Type") + xlab("") + theme(legend.position = "none",
                                                                        axis.text.x = element_text(face="bold",size =10, colour = 'black'), 
                                                                        axis.title.y = element_text(face="bold", size =10),
                                                                        axis.text.y = element_text(size =10, colour = 'black'))

regPercents_BZNY <- data.frame(comparison = c("Maternal","Maternal", "Maternal", "Maternal", "Genetic", "Genetic","Genetic","Genetic"),
                               reg_type = c("cis only", "trans only", "cis + trans", "cis x trans","cis only", "trans only", "cis + trans", "cis x trans"),
                               reg_percent = c(32.4,28.4,29.7,9.5,   52.8,17.9,27.2,2.1))

ggplot(regPercents_BZNY, aes(fill=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")), y=reg_percent, x=comparison, 
                             alluvium=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")), 
                             stratum=factor(reg_type, levels = c("cis only", "trans only", "cis + trans", "cis x trans")))) + 
  geom_flow(alpha=0.6,width = 0.6)+ geom_stratum(color = "transparent",alpha=0.95, width = 0.6)+
  geom_text(stat = "stratum",aes(label = percent(after_stat(prop), accuracy = .1)), size = 3) + scale_fill_manual("Regulation",values = c("#f95738", "#9966cd","#ffae42", "#03c04a")) +
  scale_x_discrete(limits=c("Genetic","Maternal")) +theme_classic() + ylab("Proportion of Regulatory Type") + xlab("") + theme(legend.position = "none",
                                                                                                                             axis.text.x = element_text(face="bold",size =10, colour = 'black'), 
                                                                                                                             axis.title.y = element_text(face="bold", size =10),
                                                                                                                             axis.text.y = element_text(size =10, colour = 'black'))


