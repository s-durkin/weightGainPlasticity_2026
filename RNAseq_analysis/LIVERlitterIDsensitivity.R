#LIVER litter ID sensitivity testing


############################################################## #
#### Litter-level DESeq2 sensitivity analysis ####
#### 50 iterations: one randomly selected pup per litter ####
############################################################## #

library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("~/Documents/Berkeley/Research/Maternal Effects/Gene Expression/R Scripts/")
load("liver_DE_comps.RData")

# Original DESeq2 results
# pureline_NYvsBZ = NYNY vs BZBZ
# NY_pureVcf = NYNY vs NYNY-CF
# BZ_pureVcf = BZBZ-CF vs BZBZ

# Define original significant DEGs
sigDEGsFC05_NNvBB <- pureline_NYvsBZ_liver %>%filter(padj <= 0.05 & abs(log2FoldChange) > 0.5)
sigDEGsFC05_NNvNNCF <- NY_pureVcf_liver %>%filter(padj <= 0.05 & abs(log2FoldChange) > 0.5)
sigDEGsFC05_BBCFvBB <- BZ_cfVpure_liver %>%filter(padj <= 0.05 & abs(log2FoldChange) > 0.5)


############################################################## #
#### Prepare metadata ####
############################################################## #

liver_meta <- as.data.frame(colData(dds_liver))

# Make litter IDs unique across experimental groups
liver_meta$group <- interaction(liver_meta$genotype,liver_meta$maternal,drop = TRUE)

liver_meta$litter <- interaction(liver_meta$group,liver_meta$litterID,drop = TRUE)


############################################################## #
#### Run 50 iterations ####
############################################################## #

#choose one pup from each litter and re-run differential expression analysis, 50x

set.seed(123)
n_iter <- 50

# Store results from each iteration
results_NNvBB_litterTest <- vector("list", n_iter)
results_NNvNNCF_litterTest <- vector("list", n_iter)
results_BBCFvBB_litterTest <- vector("list", n_iter)

# Store selected sample IDs
selected_samples <- vector("list", n_iter)

run_times <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  
  cat("Running iteration", i, "of", n_iter, "\n")
  start_time <- Sys.time()
  
  # Randomly choose one pup from each litter
  selected <- unlist(lapply(split(rownames(liver_meta), liver_meta$litter),
                            function(x) sample(x, 1)),
                     use.names = FALSE)
  
  selected_samples[[i]] <- selected
  
  # Subset original DESeq2 object
  dds_liver_sub <- dds_liver[, selected]
  
  # Remove unused factor levels
  dds_liver_sub$genotype <- droplevels(dds_liver_sub$genotype)
  dds_liver_sub$maternal <- droplevels(dds_liver_sub$maternal)
  
  # Original design
  design(dds_liver_sub) <- ~ genotype + maternal + genotype:maternal
  
  # Get model matrix
  mod_mat_liver <- model.matrix(
    design(dds_liver_sub),
    colData(dds_liver_sub))
  
  # Define coefficient vectors
  NY_pureline <- colMeans(mod_mat_liver[
    dds_liver_sub$genotype == "NY" &
      dds_liver_sub$maternal == "NY", , drop = FALSE])
  
  BZ_pureline <- colMeans(mod_mat_liver[
    dds_liver_sub$genotype == "BZ" &
      dds_liver_sub$maternal == "BZ", , drop = FALSE])
  
  NY_cf <- colMeans(mod_mat_liver[
    dds_liver_sub$genotype == "NY" &
      dds_liver_sub$maternal == "BZ", , drop = FALSE])
  
  BZ_cf <- colMeans(mod_mat_liver[
    dds_liver_sub$genotype == "BZ" &
      dds_liver_sub$maternal == "NY", , drop = FALSE])
  
  # Run DESeq2
  deseq_liver_sub <- DESeq(dds_liver_sub, test = "Wald",quiet = TRUE)
  
  # NYNY vs BZBZ
  res_pureline_NYvsBZ <- results(deseq_liver_sub,
                                 contrast = NY_pureline - BZ_pureline,pAdjustMethod = "BH")
  
  # NYNY vs NYNY-CF
  res_NY_pureVcf <- results(deseq_liver_sub,
                            contrast = NY_pureline - NY_cf,pAdjustMethod = "BH")
  
  # BZBZ-CF vs BZBZ
  res_BZ_cfVpure <- results(deseq_liver_sub,
                            contrast = BZ_cf - BZ_pureline ,pAdjustMethod = "BH")
  
  # Save results
  results_NNvBB_litterTest[[i]] <- as.data.frame(res_pureline_NYvsBZ)
  results_NNvNNCF_litterTest[[i]] <- as.data.frame(res_NY_pureVcf)
  results_BBCFvBB_litterTest[[i]] <- as.data.frame(res_BZ_cfVpure)
  
  run_times[i] <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  
  cat("Completed in", round(run_times[i], 1), "seconds\n")
}



############################################################## #
#### Compare original DEGs with subsampled results ####
############################################################## #

# Extract log2FC for the ORIGINAL significant DEG sets
genes_NNvBB <- rownames(sigDEGsFC05_NNvBB)
genes_NNvNNCF <- rownames(sigDEGsFC05_NNvNNCF)
genes_BBCFvBB <- rownames(sigDEGsFC05_BBCFvBB)

# Original fold changes
original_LFC_NNvBB <-pureline_NYvsBZ_liver[genes_NNvBB, "log2FoldChange"]
original_LFC_NNvNNCF <-NY_pureVcf_liver[genes_NNvNNCF, "log2FoldChange"]
original_LFC_BBCFvBB <-BZ_cfVpure_liver[genes_BBCFvBB, "log2FoldChange"]


# Matrix: genes x iterations
lfc_matrix_NNvBB <- sapply(results_NNvBB_litterTest,
                           function(x) x[genes_NNvBB, "log2FoldChange"])

lfc_matrix_NNvNNCF <- sapply(results_NNvNNCF_litterTest,
                             function(x) x[genes_NNvNNCF, "log2FoldChange"])

lfc_matrix_BBCFvBB <- sapply(results_BBCFvBB_litterTest,
                             function(x) x[genes_BBCFvBB, "log2FoldChange"])

rownames(lfc_matrix_NNvBB) <- genes_NNvBB
rownames(lfc_matrix_NNvNNCF) <- genes_NNvNNCF
rownames(lfc_matrix_BBCFvBB) <- genes_BBCFvBB


############################################################## #
#### Correlation for each iteration ####
############################################################## #

cor_NNvBB <- apply(lfc_matrix_NNvBB, 2,function(x) cor(original_LFC_NNvBB, x,use = "complete.obs"))
cor_NNvNNCF <- apply(lfc_matrix_NNvNNCF, 2,function(x) cor(original_LFC_NNvNNCF, x,use = "complete.obs"))
cor_BBCFvBB <- apply(lfc_matrix_BBCFvBB, 2,function(x) cor(original_LFC_BBCFvBB, x,use = "complete.obs"))

#info on the individual correlations - lowest in the iterations, median, mean, maximum
summary(cor_NNvBB) #0.9909 - 0.9940, avg. 0.9927

summary(cor_NNvNNCF) #0.9714 - 0.9870, avg. 0.9816

summary(cor_BBCFvBB) #0.9481 - 0.9791, avg. 0.9668



############################################################## #
#### Average log2FC across 50 iterations ####
############################################################## #

mean_LFC_NNvBB <- rowMeans(lfc_matrix_NNvBB, na.rm = TRUE)
mean_LFC_NNvNNCF <- rowMeans(lfc_matrix_NNvNNCF, na.rm = TRUE)
mean_LFC_BBCFvBB <- rowMeans(lfc_matrix_BBCFvBB, na.rm = TRUE)

compare_NNvBB <- data.frame(gene = genes_NNvBB,original_LFC = original_LFC_NNvBB,mean_subsample_LFC = mean_LFC_NNvBB)
compare_NNvNNCF <- data.frame(gene = genes_NNvNNCF,original_LFC = original_LFC_NNvNNCF,mean_subsample_LFC = mean_LFC_NNvNNCF)
compare_BBCFvBB <- data.frame(gene = genes_BBCFvBB,original_LFC = original_LFC_BBCFvBB,mean_subsample_LFC = mean_LFC_BBCFvBB)



############################################################## #
#### Plot correlations ####
############################################################## #

plot_litter_cor <- function(df, title) {
  
  r <- cor(df$original_LFC,df$mean_subsample_LFC,use = "complete.obs")
  
  ggplot(df, aes(x = original_LFC,y = mean_subsample_LFC)) +
    geom_point(alpha = 0.35, size = 1.5) +geom_abline(slope = 1, intercept = 0,linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE,color = "red") +
    annotate("text", x = -Inf, y = Inf,label = paste0("Pearson r = ", round(r, 3)),hjust = -0.1, vjust = 1.5) +
    theme_classic() +labs(title = title,x = "Original DESeq2 log2FC",y = "Mean log2FC across 50 litter subsamples")
}

plot_litter_cor(compare_NNvBB, "SARA vs MANA liver")
plot_litter_cor(compare_NNvNNCF, "SARA vs SARA-CF liver") + coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4))
plot_litter_cor(compare_BBCFvBB, "MANA-CF vs MANA liver") + coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))






